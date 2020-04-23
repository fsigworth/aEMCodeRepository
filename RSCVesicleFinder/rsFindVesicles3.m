function [mi, state]=rsFindVesicles3(m,mi,pars,findInMask)
% Auto vesicle finder for Vesicle_finding_GUI
% Calls:
% [mi t]=meFindVesicles3(m, mi, pars, findInMask)          --initialize vesicle finder
% meFindVesicles('end');      --deallocates the persistent variables.
% mi = meFindVesicles3('next',maxN,minmaxThresh);  --find up to maxN vesicles.
% mi = meFindVesicles3(m, mi) -- new image, old parameters
% 
% pars is a struct with fields
%   ds0 : scale up to global coords
%   ds0Shift : shift to subtract from global coords
%   rPars : min, max, step for radii
% t is a persistent structure that contains information for the
% The returned mi.vesicle.ok field for each found vesicle
% has the value [1 f 0 0] where f=1 for a
% vesicle entirely in range, f=0 if it falls a bit below the amplitude
% threshold or above the max radius.

persistent t;

useSpheres=1;
% blankRadiusFactor=.7;
borderFactor=1.25;
maxFracMasked=.6;   % tries to fit vesicles up to this overlap with mask.
maxGoodFracMasked=.5;
extendedFitRadius=1.6;
reducedFitRadius=0.8;
radiusThresh=.95;  % only accept vesicles larger than this fraction of rmin.
%  bigRThresh=.5; % asymptotic rDecay value, now given by "ice amp", amps(3).

% targetPixA=10;
displayOn=0;
rExponent=.5;  % radius-weighting

if isnumeric(m)  % we are initializing the finder.
    t.mi=mi;
    ns=size(m); % use the display image size
    n0=round(mi.imageSize/pars.ds0); % expected size of the stored mask
    msk=Crop(meGetMask(mi,n0),ns); % convert msk to logical

    m=m-median(m(logical(msk(:))));  % subtract the median of all unmasked points
    m=m.*msk;
    
    % Get image and pixel sizes
%     ds0=mi.imageSize(1)/ns(1);  % downsampling factor of m
    nsx=NextNiceNumber(borderFactor*ns);  % expanded number of pixels
    t.borderShift=(nsx-ns)/2;  % displacement of the center of padded image from unpadded center
    t.msx=Crop(m,nsx);  % pad the image
    t.ms=m;
    t.ns=ns;
    t.ds=pars.ds0;
    t.dsShift=pars.ds0Shift;
    pixA=t.ds*mi.pixA;  %pixA in the image ms
    t.pixA=pixA;
    
    sHP=.002*pixA;  % highpass frequency 500 A
    varLP=.004*pixA;  % lowpass for image 250 A

    vm=meDownsampleVesicleModel(mi.vesicleModel,t.ds)*pixA;
    %   note that the density is scaled up to match the scattering per voxel.
    t.mbnThickness=(sum(vm>max(vm)/4));  % approx. membrane thickness in pixels.
    
%     Make a Gaussian HP filter
    highPass=exp(-sHP^2./(max(RadiusNorm(nsx),1e-6).^2));  % zero-center highpass function
    t.fms=fftn(t.msx).*ifftshift(highPass);  % Fourier-transformed, filtered image.
    
    if displayOn
        subplot(2,3,1);
        imags(t.ms);
    end;
    
    % Get the effective CTF from the merging.
%     disp('Finder start: getting effCTF');

    H=ifftshift(highPass.*meGetEffectiveCTF(mi,nsx,t.ds));
    
    % Get the mask, pad it to the search size
    t.mask=Crop(msk,nsx);
%     t.mask=Crop(meGetMask(mi,ns),nsx);
    t.findInMask=findInMask;
    
    %%
    rPars=pars.rPars;
    if ~(isfield(t,'rPars') && all(t.rPars==rPars) ...
            && nsx(1) == size(t.frefs,1))  % need to compute references
        disp('Constructing references');
        t.rPars=rPars;
        t.rMin=rPars(1)/pixA;  % minimum liposome radius in pixels
        t.rMax=rPars(2)/pixA;  % maximum liposome radius
        t.fitmin=t.rMin*reducedFitRadius;
        t.fitmax=t.rMax*extendedFitRadius;
        t.rstep=rPars(3)/pixA; % radius step
        t.nrsteps=round((t.fitmax-t.fitmin)/t.rstep+1);
        
        t.frefs=single(zeros([nsx,t.nrsteps]));

        t.fsph=single(zeros([nsx,t.nrsteps]));
        t.spowers=zeros(t.nrsteps,1);
        
        t.localSD=zeros([nsx,t.nrsteps],'single');
        t.rads=zeros(t.nrsteps,1);
        t.powers=zeros(t.nrsteps,1);
        t.rrefs=single(zeros([nsx t.nrsteps]));
        t.orefs=single(zeros([nsx t.nrsteps]));
        t.opowers=zeros(t.nrsteps,1);
        
        if displayOn            
            subplot(2,3,2);
            nd=NextNiceNumber(3*t.rMax);  % size of box to display
            nd=min(nd,nsx);
        end;
        vesOrigin=floor(nsx/2+1);   % shift the origin of the references
%         so they'll correspond to the center of the padded image.  Zero
%         lag will then be (1,1), pointing to the correct location int the
%         padded image.
        
        fmVar=fftn(GaussFilt(t.msx,varLP).^2);
        for i=1:t.nrsteps
%             Make vesicle model
            r=t.fitmin+(i-1)*t.rstep;  % radius in pixels
            t.rads(i)=r;  % radius in A
            v = VesicleFromModel(nsx,r,vm,vesOrigin);
                 
%             Compute the local SD, ns points
            rd=ifftshift(fuzzymask(nsx,2,r+1.5*t.mbnThickness,t.mbnThickness/2)); %%% arbitrary
            fd=fftn(rd); % disc and ft for variance calc
            lsd=sqrt(abs(real(ifftn(fmVar.*fd))/(rd(:)'*rd(:))));
            t.localSD(:,:,i)=lsd;

            fv=-H.*fftn(ifftshift(v));    % FT of vesicle at origin.
            t.frefs(:,:,i)=fv; % FT of filtered, CTF vesicle
            t.orefs(:,:,i)=v;  % original vesicle
            t.opowers(i)=v(:)'*v(:); % original vesicle variance.
            
            rv=fftshift(real(ifftn(fv)));  % CTF-filtered real image of vesicle, centered
            t.powers(i)=rv(:)'*rv(:); % Filtered vesicle variance
            t.rrefs(:,:,i)=rv;        % Filtered vesicle

            if useSpheres  % make inner sphere density too
                sphScale=2*pixA;  % assume 2 volt inner potential
                sph=sphScale*SphereDensity(nsx,r-t.mbnThickness/2,vesOrigin);
                fs=-H.*fftn(ifftshift(sph));
                t.fsph(:,:,i)=fs;     % Corresponds to t.frefs
                rsph=fftshift(real(ifftn(fs)));
%  imags(rsph);
%  drawnow;
                t.spowers(i)=rsph(:)'*rsph(:); % corresponds to t.powers
            end;
            if displayOn
                %                 t.rrefs(:,:,i)=rv;
                imags(rv);
                title(i);
                drawnow;
            end;
        end;
    end;
% %     medSD=median(t.localSD(:)); % ??median over all radii??
%     Make localSD bigger.
    t.localSD=Crop(Crop(t.localSD,ns,1),nsx,1,1); % SD=1 outside the original image
    %%
    % Compute all the cross-correlations, using a weight of r^rExponent to equalize
    % the CC peak values
    disp('Computing cross-correlations');
    t.ccs=single(zeros([nsx t.nrsteps]));
    t.nccs=single(zeros([nsx t.nrsteps]));
    t.spccs=zeros([nsx t.nrsteps],'single');
    % cc values are scaled by 1/r^rExponent to make small vesicles as easily found
%     rExpNorm=t.rMin^rExponent;
    for i=1:t.nrsteps
        rScale=(t.rads(i)/t.rMin)^(rExponent); % =1 at rMin, increases with rExponent
%       t.ccs is scaled by 1/t.power so its value is the amplitude.
        t.ccs(:,:,i)=rScale/t.powers(i)*...
            real(ifftn(t.fms.*conj(t.frefs(:,:,i))));
% %         t.nccs(:,:,i)=t.ccs(:,:,i)./t.localSD(:,:,i)*t.rads(i)^-rExponent;
        t.spccs(:,:,i)=real(ifftn(t.fms.*conj(t.fsph(:,:,i))))/t.spowers(i); % sphere CC
    end;
    [t.ccsmx, t.ccsmi]=max(t.ccs,[],3);
    t.radScalings=((t.rads/t.rMin).^-rExponent);  % undo scaling of references
    t.ccsmxScaled=t.ccsmx.*t.radScalings(t.ccsmi); % true scaling
%     We'll determine the amplitude and final radius assignments made from
%     rescaled ncc functions.
    ccsScaled=zeros([nsx t.nrsteps],'single');
    for k=1:t.nrsteps
        ccsScaled(:,:,k)=t.ccs(:,:,k)*t.radScalings(k);
    end;
    [t.ccsmxScaled,t.ccsmiScaled]=max(ccsScaled,[],3); % we'll use these for manual picking.
    t.sphMxCC=max(t.spccs,[],3);
    
    %%
    disp('Finding vesicles');
    
    if displayOn
        subplot(2,3,2);
        imags(t.ccsmx.*t.mask);
        drawnow;
    end;
    if findInMask
        t.ccsmx2=t.ccsmx;
    else
        t.ccsmx2=t.ccsmx.*t.mask;
    end;
    
    if displayOn
        t.msub=t.ms;
    end;
    t.umodel=zeros(ns);
    
    %    % Check to see how many entries are already there.
    t.nfound=numel(t.mi.vesicle.x);
    
%%-----------------------------------Finding loop--------------------------------    
else % m is a string, 'next' or 'end'
    switch lower(m)
        case 'next'
            maxN=mi;        % pick up the alternate arguments.
            thresh=pars;
            bigRThresh=pars(3); % allowed decay for large vesicles
            thresh(3)=0; % ice amp is constrained to this.
            thresh(4)=thresh(1) * 0.5;  % arbitrary factor!!!
            nsx=size(t.ccsmx2);
            ctrx=floor(nsx/2+1); % distance to origin
% %             medSD=median(t.localSD(:));

            nf0=t.nfound;
            while t.nfound<nf0+maxN
                [t.globalmax, jx, jy]=max2d(t.ccsmx2);
                if t.globalmax<1e-6 % reasonable threshold should be >1e-4
                    break
                end;
                % Get the corrected peak value and interpolated radius index
                [ampi, refii]=max1di(squeeze(t.ccs(jx,jy,:)).*t.radScalings);
                jz=max(1,min(round(refii),t.nrsteps)); % nearest index for correct radius vesicle
                refri=t.fitmin+(refii-1)*t.rstep;  % interpolated radius in pixels
                rDecay=bigRThresh+(1-bigRThresh)/(1+(refri*t.pixA/200)^2); % help for big vesicles.

%                     Blank the cc peak and check the fraction masked
%                 refi=max(1,min(round(refii),t.nrsteps)); % closest model radius
                support=single(circshift(t.orefs(:,:,jz),round([jx jy])-ctrx) ...
                                                /sum(t.mi.vesicleModel) > .001);
                t.ccsmx2=t.ccsmx2.*(~support); % blank the running cc map
%  imags(t.ccsmx2); drawnow;
                fracMasked=support(:)'*single((~t.mask(:)))/sum(support(:));
                if max(ampi,t.globalmax)>thresh(1)*(1-maxFracMasked)*bigRThresh % quick check for any possible find
                    spcc=t.spccs(jx,jy,jz); % Check the sphere correlation too.
% %                 nccv=t.nccs(jx,jy,jz)*medSD;
% %                 blank=fuzzymask(nsx,2,blankRadiusFactor*refri+t.mbnThickness/2,t.mbnThickness,[jx jy]);

                %                 Test for a found vesicle.
                    if fracMasked<maxFracMasked ...
                        && ampi>thresh(1)*(1-fracMasked)*rDecay...
                         && spcc<thresh(3)...
                        && (refri>=radiusThresh*t.rMin)
%                   disp([num2str(t.nfound+1) '  ' num2str(ampi) '  spcc:' num2str(spcc) '  r:' num2str(refri)...
%                       '  fracMasked:' num2str(fracMasked) '  msk:' num2str(t.mask(jx,jy))...
%                       '  rDecay:' rDecay]);
                        % We found something. Now set a flag
                        flag=single((fracMasked<maxGoodFracMasked)...
                            && (refri>=radiusThresh*t.rMin)...
                            && (refri<=t.rMax)...
                            && t.mask(jx,jy));
%                         && nccv>thresh(4);

                        t.nfound=t.nfound+1; %%%%% increment nfound
                        t.mi.vesicle.r(t.nfound,1)=refri*t.ds;
%                         We allow vesicle x and y to be out of bounds.
                        t.mi.vesicle.x(t.nfound,1)=(jx-1-t.borderShift(1))*t.ds+1-t.dsShift(1);
                        t.mi.vesicle.y(t.nfound,1)=(jy-1-t.borderShift(2))*t.ds+1-t.dsShift(2);
                        t.mi.vesicle.s(t.nfound,1)=ampi/(1-fracMasked);
                        t.mi.vesicle.ok(t.nfound,1:4)=[1 flag 0 0];  % flag indicates a vesicle in bounds.
                        vref=ampi*circshift(t.rrefs(:,:,jz),round([jx jy]-nsx/2-1));
                        t.umodel=t.umodel+Crop(vref,t.ns);  % approximate model
                    end;
                else
                    break;
                end;
            end;
            if displayOn
                subplot(2,3,4); imags(t.msub);
                title(t.nfound);
                subplot(2,3,2); imags(t.umodel);
                subplot(2,3,6); imags(t.ccsmx2);
                
                subplot(2,3,3);
                plot(t.mi.vesicle.r(:,1)*t.mi.pixA,t.mi.vesicle.s(:,1),'k.');
                xlabel('Vesicle radius, �');
                ylabel('Image amplitude');
                drawnow;
            end;
            % %%
            % % Create new vesicle model
            % model=zeros(ns,ns);
            % for i=1:nfound
            %     model=model+amps(i)*VesicleDensity(ns,radii(i),mbnThickness,coords(:,i)+1);
            % end;
            
        case 'end'
            t=[];  % clear the whole struct
    end;  % switch
end;  % m numeric
t.mi.vesicle.refined=0;
mi=t.mi;
state=t;


%%
%
%     disp('Final subtraction of vesicles');
%     vm=meMakeModelVesicles(mi,n);
%     mv=m-vm;
%     %%
%     figure(2);
%     imacs(GaussFilt(mv,.2));
%
%     figure(1);
%     subplot(2,3,1);
%     imacs(BinImage(m,4));
%     subplot(2,3,2);
%     imacs(vm);
%     subplot(2,3,3);
%     plot(mi.vesicle.r*mi.pixA,mi.vesicle.s,'k.','markersize',10);
%
%     subplot(2,3,4);
%     imacs(BinImage(mv,4));
%     title('Subtracted');
%     subplot(2,3,5);
%     hist(mi.vesicle.s,50);
%     xlabel('Vesicle amplitude s');
%     drawnow;
%
%     mi.basePath=ParsePath(inPath);  % make it the local path
%     bname=[mi.basePath mi.procPath mi.baseFilename];
%     WriteMRC(mv,ds0*mi.pixA,[bname 'mv.mrc']);
%     jname=[mi.basePath mi.procPath 'jpeg/' mi.baseFilename];
%     imwrite(uint8(imscale(rot90(mv),256,1e-3)),[jname 'mv.jpg']);
%
