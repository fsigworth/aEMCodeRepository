function mi1=rsRefineVesicleFitsSub(mi,m,fittingMode)
% modes are:  1 - fit overall amplitude
%             2 - fit amplitude and shift only, don't change radius
%             3 - fit amplitude, shift and radius
if nargin<3
    fittingMode=3;
end;
maxVesiclesToFit=inf;
useOkField=1;      % refine every vesicle for which ok is true.
doDownsampling=1;  % Downsample for speed
disA=800;          % size of displayed/fitted window in angstroms
sAcceptSd=0.6;       % number of SDs to allow in ok vesicles
% maxRadius=400;     % maximum radius in Å for ok to be true.
% minRadius=100;
% fileName='';
displayFreq=10;
%         sp1=RadialPowerSpectrum(m);  % This is very slow.
%         nsp=numel(sp1);
%         hfVar0=mean(sp1(round(0.3*nsp):round(0.45*nsp)))
%       faster code:
        n=size(m);
        annulus=fuzzymask(n,2,0.225*n,.05*n)-fuzzymask(n,2,0.15*n,.05*n);
        spc=annulus.*fftshift(abs(fftn(m)).^2)/(n(1)*n(2));
        hfVar0=sum(spc(:))/sum(annulus(:));
        

        % Get image and pixel sizes
        n=size(m,1);
        ds0=mi.imageSize(1)/n;  % downsampling factor of m
        pixA0=mi.pixA*ds0;    % pixel size of m
        if doDownsampling
            % downsample to about 10A per pixel, yielding the image ms
            targetPixA=12;  % maximum pixel size
            ns=NextNiceNumber(n*pixA0/targetPixA,5,4);  % multiple of 4, max factor 5.
            if ns<n
                disp(['Downsampling to ' num2str(ns) ' pixels.']);
                ms=Downsample(m,ns);
            else
                ns=n;
                ms=m;
            end;
            ds=ds0*n/ns;  % downsampling factor of ms relative to original images.
            pixA=ds*mi.pixA;  % pixA in the image ms.
        else  % use the original merged image scale
            ds=ds0;
            pixA=pixA0;
            ns=n;
            ms=m;
        end;
        hfVar=hfVar0*(ds0/ds)^2;  % hf spectral density after downsampling
        ndis=NextNiceNumber(disA/pixA);  % size of display/fitting image
        
        %%  Get the original subtraction, and modify the amplitudes if necessary.
        vs=meMakeModelVesicles(mi,ns);  % old subtraction
        msAmpScale=(vs(:)'*ms(:))/(vs(:)'*vs(:))
        if msAmpScale>1e-3 % don't allow ridiculous values
            vs=vs*msAmpScale;
            mi.vesicle.s=mi.vesicle.s*msAmpScale;
        end;
        msub=ms-vs;
        mmask=meGetMask(mi,ns);
%         Whiten the subtracted image
        pwH=meGetNoiseWhiteningFilter(mi,ns);
        msubf=mmask.*real(ifftn(fftn(msub).*ifftshift(pwH)));
        figure(2);
        SetGrayscale;
        imacs(msubf);
        figure(1);
        SetGrayscale;
        subplot(2,3,1);
        imacs(msubf);
        title(['Old subtraction, scaling ' num2str(msAmpScale)]);

        mi1=mi;  % we'll load the new coordinates into mi1

        
        %%
        if fittingMode>1  % tune up the x, y and s for each vesicle.
                %
                if ~isfield(mi.vesicle,'ok')
                    mi1.vesicle.ok=ones(size(mi1.vesicle.x));
                end;
                %%  Actual fitting is done here
                figure(1)
                disp('Fine fitting translation and amplitude.');
                drawnow;
                nVesicles=sum(mi1.vesicle.ok>0);
                nVesicles=min(nVesicles,maxVesiclesToFit);
                nVesicles
                effCTF=ifftshift(meGetEffectiveCTF(mi,ndis,ds));
                pwFilter=ifftshift(meGetNoiseWhiteningFilter(mi,ndis,ds));
                figure(1);
                for ind=1:nVesicles
                    displayOn=mod(ind,displayFreq)==0;  % Display only some images
                    ok=mi1.vesicle.ok(ind)>0;
                    if ~useOkField || ok
                        [mi1 diffIm vesFit]=rsQuickFitVesicle(msubf,mmask,mi1,mi,...
                            ind,effCTF.*pwFilter,hfVar,fittingMode,displayOn);
                        
                        subplot(2,2,2); title(mi.baseFilename,'interpreter','none');
                    end;
                end;
                maxR=maxRadius/mi1.pixA;
                minR=minRadius/mi1.pixA;
                medianS=median(mi1.vesicle.s);
                mi1.vesicle.ok=(mi1.vesicle.s>medianS/(1+sAcceptSd)) & (mi1.vesicle.s<(1+sAcceptSd)*medianS)...
                    & (mi1.vesicle.r < maxR) & (mi1.vesicle.r > minR);
                numberOK=sum(mi1.vesicle.ok)
                mi1.vesicle.refined=1;
         end; % switch

