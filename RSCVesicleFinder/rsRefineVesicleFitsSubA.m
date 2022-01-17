function miNew=rsRefineVesicleFitsSubA(miOld,m,pars)
% Do the vesicle refinement.  miOld is the original mi file; miNew is a
% copy with model parameters, vesicle.extraPeaks and vesicle.extraSD set. m
% is the original image, described by the transform matrix pars.M1

ps=pars; % structure to pass to sub-function for small fits

minRadiusA=45;
nZeros=1;          % number of zeros used in the merged CTF
tinySValue=1e-4;   % negligible amplitude threshold

%       Calculate variance from .30 to .45 x Nyquist; this is faster than using RadialPowerSpectrum:
ns=size(m);       % Size of the 8A image.
% pick the frequency range with an annulus in freq domain
annulus=fuzzymask(ns,2,0.225*ns,.05*ns)-fuzzymask(ns,2,0.15*ns,.05*ns);
spc=annulus.*fftshift(abs(fftn(m)).^2)/(prod(ns));
hfVar=sum(spc(:))/sum(annulus(:));
pixAs=pars.M(1,1)*miOld.pixA; % pixel size of the image we're given.
% % Handle cases where the ok field is old style
% if ~isfield(miOld.vesicle,'ok') || numel(miOld.vesicle.ok)<numel(miOld.vesicle.x) % no ok field at all
%     miOld.vesicle.ok=true(numel(miOld.vesicle.x),1);
% end;

nv=numel(miOld.vesicle.x);

% Create the new mi structure
miNew=miOld;

% Insert the new 'active fraction' field.
if ~isfield(miNew.vesicle,'af') || numel(miNew.vesicle.af)~=nv  % active fraction (= unmasked)
    miNew.vesicle.af=ones(nv,1,'single');
end;

%%  Get the original subtraction, and modify the amplitudes if necessary.
% miNew should be pretty much the same as miOld except for the fields
% vesicle.extraPeaks and vesicle.extraSD.
%

sVals=miOld.vesicle.s(:,1,1);  % Copy the basic amplitude
sVals(isnan(sVals)|sVals<tinySValue)=0;
miNew.vesicle.ok(sVals==0,:)=false;  % zero or negative assigned to be null vesicles.
miOld.vesicle.ok(sVals==0,:)=false;

finalTerms=numel(pars.rTerms); % maximum number of terms allowed.
pars.finalTerms=finalTerms;

% Make vesicle.r have the correct number of elements
miNew.vesicle.r(:,end+1:finalTerms)=0;
miNew.vesicle.r(:,finalTerms+1:end)=[];


vesList=find(miNew.vesicle.ok(:,1));
nvToFit=min(numel(vesList),pars.maxVesiclesToFit);
disp([num2str(nvToFit) ' vesicles to fit.']);

scls.M=pars.M;
scls.n=size(m);
dss=scls.M(1,1);

%   Compute the old subtraction (small image size)
nsPW=meGetNoiseWhiteningFilter(miOld,ns,dss,nZeros,pars.fHP);
nsCTF=meGetEffectiveCTF(miOld,ns,dss);
msf=real(ifftn(fftn(m).*ifftshift(nsPW)));  % High-pass filtered image

if pars.displayOn
    figure(1);
    clf;
    imags(GaussFilt(msf,.1*dss));
    title(['Prewhitened image ' miOld.baseFilename],'interpreter','none');
    drawnow;
end;

%     Make a subtraction of the entire image
if pars.doPreSubtraction
    disp('Subtracting the old vesicle fits.')
    vs=meMakeModelVesicles(miOld,scls,vesList,0,0);
    vsf=real(ifftn(fftn(vs).*ifftshift(nsPW.*nsCTF)));
    if pars.displayOn
        imags(GaussFilt(m-vsf,.1*pixAs));
        title(['Preliminary subtraction ' miOld.baseFilename],'interpreter','none');
        drawnow;
        figure(2);
    end;
else
    vsf=0;
end;

% Get the CTF information for the 1-vesicle fitting regions
nds=NextNiceNumber(pars.disA/pixAs);  % size of display/fitting image
ndsPW=meGetNoiseWhiteningFilter(miOld,nds,dss,nZeros,pars.fHP);
ndsCTF=meGetEffectiveCTF(miOld,nds,dss);

%         Get the masks
layers=1:min(pars.maxMaskLayers,numel(miOld.mask));
% msmask=Crop(meGetMask(miOld,round(miOld.imageSize/dss),layers),ns);
msmask=meGetMask(miOld,ns,layers);

if pars.fitRadius % we're doing nonlinear fit


    %%  Actual radius fitting is done here
    if pars.listFits
        disp('  ind 1000s  ---r (A)----   pick   ok   nTerms ------- 100s/s(1) -------------');
        %        1  2.779  205   0  20     2  1 1 1 1    4    24.80   23.71    0.00    0.00   0
    end;

    miNew.vesicle.ok(:,3)=miNew.vesicle.ok(:,1);  % we'll mark unfitted vesicles here.

    ps=struct;
    ps.displayOn=pars.displayOn;

    for j=1:nvToFit  % Loop over vesicles
        ind=vesList(j);

        %             Look up the number of radius terms to use.
        %             rTerms is something like [150 200 250 300 350 400 inf];
        %             r<150� gets 1 term; <200 gets 2 terms; etc.
        % Set the number of amplitude terms.

        %-------------------Basic fit------------------

        % Look up the number of radius terms to fit

        ps.finalTerms=find(miNew.vesicle.r(ind,1)<pars.rTerms/miNew.pixA,1);


        %%%% rConstraints set here.
        ps.rConstraints=ones(1,ps.finalTerms);
        ps.rConstraints(2:ps.finalTerms)=0.4./((2:ps.finalTerms).^2)';
        ps.preservedTerms=pars.preservedTerms;
        ps.nIters=pars.nIters;
        ps.fHP=pars.fHP;
        ps.hfVar=hfVar;
        ps.extraRound=pars.extraRound;
        if pars.doPreSubtraction
            % First, compute the one vesicle in question
            vs1=meMakeModelVesicles(miOld,scls,ind,0,0); % no ctf or prewhitening
            vs1f=real(ifftn(fftn(vs1).*ifftshift(nsPW.*nsCTF)));  % filtered model
        else
            vs1f=0;
        end;
        ps.M=scls.M;
        %         Repeated fits with perturbed initial radius
        ndr=numel(pars.radiusStepsA);
        miTemps=cell(ndr,1);
        errs=zeros(ndr,1);
        resImgs=zeros(nds,nds,ndr,'single');
        rSteps=pars.radiusStepsA/miNew.pixA;

        for jr=1:ndr
            %             Perturb the radius
            miInput=miNew;
            newR1=max(miInput.vesicle.r(ind,1)+rSteps(jr),minRadiusA/miInput.pixA);
            miInput.vesicle.r(ind,1)=newR1; % Replace the constant radius

            %      --------------nonlinear fitting---------------

            [miTemps{jr}, fitIms, vesFits]=rsRefineVesRadius(msf-vsf,vs1f,msmask,miInput,ind,ndsCTF.*ndsPW,ps);
            %  miTemps{jr},fitIms,vesFits]=rsQuickFitVesicle2(msf-vsf,vs1f,msmask,miInput,...
            %                 ind,ndsCTF.*ndsPW,ps,displayOn);

            errs(jr)=miTemps{jr}.vesicle.err(ind);
            %             if pars.displayOn
            %                 resImgs(:,:,jr)=fitIms-vesFits;
            %                 subplot(2,2,4);
            %                 imags(resImgs(:,:,jr));
            %                 title(num2str([jr pars.finalTerms]));
            %                 drawnow;
            %             end;
        end;
        %         Find the best fit
        [~,jr]=min(errs);
        [~,jr0]=min(abs(pars.radiusStepsA));

        miNew=miTemps{jr}; % we finally assign miNew here.

        if pars.displayOn && ndr>1 % show the various fit results
            for k=1:ndr
                subplot(6,4,k+16);
                imags(resImgs(:,:,k));
                if k==jr
                    str='***';
                else
                    str='';
                end;
                title([num2str(errs(k)) str]);
            end;
        end;
        rvals=zeros(1,3);
        str=sprintf(['%4d %6.3f %4.1d %4.1d %4.1d  %2d %2d%2d%2d%2d  %2d  '],...
            ind, 1000*miNew.vesicle.s(ind,1),...
            rvals(1),rvals(2),rvals(3),...
            jr-jr0,...
            miNew.vesicle.ok(ind,:), sum(miNew.vesicle.r(ind,:)~=0),0);
        disp(str);
    end; % for j over nv

    %% -------------Amplitude fitting done here ---------
else % ~pars.fitRadius
    ps.M=scls.M;
    ps.fHP=pars.fHP;
    ps.peakPositionA=pars.peakPositionA;
    ps.peakSigmaA=pars.peakSigmaA;

    miNew=miOld;
    miNew.vesicle.s=miOld.vesicle.s(:,1,1); % truncate all the amplitudes

    for j=1:nvToFit
        ind=vesList(j);
        %             Look up the number of amplitude terms to use.
        ps.finalTerms=find(miNew.vesicle.r(ind,1)<pars.aTerms/miNew.pixA,1);
        ps.fHP=pars.fHP;
        if pars.doPreSubtraction
            % Compute the one vesicle in question
            vs1=meMakeModelVesicles(miOld,scls,ind,0,0); % no ctf or prewhitening
            vs1f=real(ifftn(fftn(vs1).*ifftshift(nsPW.*nsCTF)));  % filtered model
        else
            vs1f=0;
        end;
        miNew=rsRefineVesAmplitude(msf-vsf,vs1f,msmask,miNew,ind,ndsCTF.*ndsPW,ps);
    end; % for j over nv

%         Disallow negative amplitudes, and nan values
% for i=1:nv
%     miNew.vesicle.s(isnan(miNew.vesicle.s))=0;
%     if any(isnan(miNew.vesicle.s(i,:)));
%         miNew.vesicle.s(i,:)=0;
%         miNew.vesicle.ok(i,:)=0;
%     end;
% end;



end;
% numberRefined=sum(miNew.vesicle.ok(:,3))
% numberGood=sum(all(miNew.vesicle.ok(:,1:3),2))  %    exists, in range, refined
% miNew.vesicle.refined=1;

return



%     %
%     %
%     %
%     %
%     %     %   ----------------- Linear fit only --------------
%     % else % do a linear fit of the vesicle
%     %     nPeaks=numel(miNew.vesicle.extraPeaks);
%     %
%     %     jr=1;
%     %     jr0=1;
%     %     if pars.doPreSubtraction
%     %         v1=meMakeModelVesicles(miOld,sclb,ind,0,0);
%     %         v1f=real(ifftn(fftn(v1).*ifftshift(nbPW.*nbCTF)));  % filtered model
%     %     else
%     %         v1=0;
%     %         v1f=0;
%     %     end;
%     %     pb=ps;
%     %     pb.nTerms=[0 ps.nTerms(end,end)];
%     %     pb.M=pars.M4;
%     %     %                 Do a fit of only the amplitude terms
%     %     [miNew,fitIm,vesFit]=rsQuickFitVesicle2(mbf-vfb,v1f,mbmask,miNew,...
%     %         ind,ndbCTF.*ndbPW,pb,pars.displayOn & ~doFitRadius);  % no display
%     %
%     %     if pars.displayOn  % update the subtracted image
%     %         subplot(2,2,4);
%     %         imags(GaussFilt(fitIm-vesFit,pixAb*.05)); % 20 A filter
%     %         drawnow;
%     %     end;
%     % else
%     %     vesFit=vesFits;
%     % end;
%     %
%     % nsTerms=size(miNew.vesicle.s,2);
%     % ampString=repmat('%6.2f  ',1,nsTerms-2);
%     % if pars.listFits
%     %     rvals=round(abs(miNew.vesicle.r(ind,:)*miNew.pixA));
%     %     rvals(5)=0;
%     %     str=sprintf(['%4d %6.3f %4.1d %4.1d %4.1d  %2d %2d%2d%2d%2d  %2d  ' ampString],...
%     %         ind, 1000*miNew.vesicle.s(ind,1),...
%     %         rvals(1),rvals(2),rvals(3),...
%     %         jr-jr0,...
%     %         miNew.vesicle.ok(ind,:), sum(miNew.vesicle.r(ind,:)~=0),...  % insert s(1)
%     %         100*abs(miNew.vesicle.s(ind,3:end,1))/miNew.vesicle.s(ind,1,1));
%     %     disp(str);
%     % end;
%     % end;
%     %
%     % %         Disallow negative amplitudes, and nan values
%     % for i=1:nv
%     %     miNew.vesicle.s(isnan(miNew.vesicle.s))=0;
%     %     if any(isnan(miNew.vesicle.s(i,:)));
%     %         miNew.vesicle.s(i,:)=0;
%     %         miNew.vesicle.ok(i,:)=0;
%     %     end;
%     % end;
%     % numberRefined=sum(miNew.vesicle.ok(:,3))
%     % numberGood=sum(all(miNew.vesicle.ok(:,1:3),2))  %    exists, in range, refined
%     % miNew.vesicle.refined=1;
%     %
%     %
