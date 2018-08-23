function [mi, state]=rsFindVesicles3(m,mi,rPars,findInMask)
% Auto vesicle finder for Vesicle_finding_GUI
% Calls:
% [mi t]=meFindVesicles3(m, mi, rPars, findInMask)          --initialize vesicle finder
% mi = meFindVesicles3('next',maxN,minmaxThresh);  --find up to maxN vesicles.
% mi = meFindVesicles3(m, mi) -- new image, old parameters
% meFindVesicles('end');      --deallocates the persistent variables.
% t is a persistent structure that contains information for the

persistent t;

blankRadiusFactor=.7;
borderFactor=1.25;
maxFracMasked=.6;   % tries to fit vesicles up to this overlap with mask.
extendedFitRadius=1.5;
sHP=.03;  % highpass frequency in downsampled pix^-1
targetPixA=10;
displayOn=0;
rExponent=.2;  % radius-weighting

if isnumeric(m)  % we are initializing the finder.
    t.mi=mi;
    m=m-mean(m(:));
    
    
    % Get image and pixel sizes
    n0=size(m);
    ds0=mi.imageSize(1)/n0(1);  % downsampling factor of m
    pixA0=mi.pixA*ds0;    % pixel size of m
    ns0=n0; % use the original image size
%     ns0=NextNiceNumber(n0*pixA0/targetPixA);  % an even number, max about 10% too big.
    ns=borderFactor*ns0;  % expanded number of pixels
    [t.msx, mag]=DownsampleGeneral(m,ns,ns0/n0);  % resample the image, downsample factor ds
    t.ms=Crop(t.msx,ns0);
    borderShift=(ns-ns0)/2;
    t.ns0=ns0;
    t.ds=ds0/mag;  % total downsampling (number > 1) from original image
    pixA=t.ds*mi.pixA;  %pixA in the image ms
    
    vm=meDownsampleVesicleModel(mi.vesicleModel,t.ds)*pixA;
    
    %   note that the density is scaled up to match the scattering per voxel.
    t.mbnThickness=(sum(vm>max(vm)/2));  % approx. membrane thickness in pixels.
    
%     Make a Gaussian HP filter
    highPass=exp(-sHP^2./(max(RadiusNorm(ns),1e-6).^2));  % zero-center highpass function
    t.fms=fftn(t.msx).*ifftshift(highPass);  % Fourier-transformed, filtered image.
    
    if displayOn
        SetGrayscale;
        subplot(2,3,1);
        imacs(t.ms);
    end;
    
    % Get the effective CTF from the merging.
    H=ifftshift(highPass.*meGetEffectiveCTF(mi,ns,t.ds));
    
    % Get the mask, pad it to the search size
    t.mask=zeros(ns);
    t.mask(1:ns0(1),1:ns0(2))=meGetMask(mi,ns0);
    t.findInMask=findInMask;
    
    %%
    
    if ~(isfield(t,'rPars') && all(t.rPars==rPars) ...
            && ns(1) == size(t.frefs,1))  % need to compute references
        disp('Constructing references...');
        t.rPars=rPars;
        t.rmin=rPars(1)/pixA;  % minimum liposome radius in pixels
        t.rmax=rPars(2)/pixA;  % maximum liposome radius
        t.fitmin=t.rmin/extendedFitRadius;
        t.fitmax=t.rmax*extendedFitRadius;
        t.rstep=rPars(3)/pixA; % radius step
        t.nrsteps=round((t.fitmax-t.fitmin)/t.rstep+1);
        
        t.frefs=single(zeros([ns,t.nrsteps]));
        t.rads=zeros(t.nrsteps,1);
        t.powers=zeros(t.nrsteps,1);
        t.rrefs=single(zeros([ns t.nrsteps]));
        t.orefs=single(zeros([ns t.nrsteps]));
        t.opowers=zeros(t.nrsteps,1);
        
        if displayOn
            
            subplot(2,3,2);
            nd=NextNiceNumber(3*t.rmax);  % size of box to display
            nd=min(nd,ns);
        end;
        vesOrigin=(ns/2+1)+borderShift;   % shift the origin of the references
        for i=1:t.nrsteps
            r=t.fitmin+(i-1)*t.rstep;  % radius in pixels
            t.rads(i)=r;  % radius in A
            v = VesicleFromModel(ns,r,vm,vesOrigin);
            fv=-H.*fftn(ifftshift(v));    % FT of vesicle at origin.
            t.frefs(:,:,i)=fv;
            t.orefs(:,:,i)=v;
            t.opowers(i)=v(:)'*v(:);
            rv=fftshift(real(ifftn(fv)));  % CTF-filtered real image of vesicle, centered
            t.powers(i)=rv(:)'*rv(:);
            t.rrefs(:,:,i)=rv;
            if displayOn
                %                 t.rrefs(:,:,i)=rv;
                imacs(rv);
                title(i);
                drawnow;
            end;
        end;
    end;
    
    %%
    % Compute all the cross-correlations, using a weight of r^rExponent to equalize
    % the CC peak values
    disp('Computing cross-correlations');
    t.ccs=single(zeros([ns t.nrsteps]));
    % cc values are scaled by 1/r^rExponent to make small vesicles as easily found
    for i=1:t.nrsteps
        t.ccs(:,:,i)=t.rads(i).^(rExponent)/t.powers(i)*...
            real(ifftn(t.fms.*conj(t.frefs(:,:,i))));
    end;
    % ccs=ccs-min(ccs(:));  % min value is zero.
    [t.ccmx, t.ccmi]=max(t.ccs,[],3);
    t.ccmxScaled=t.ccmx.*(t.rads(t.ccmi).^-rExponent);  % undo scaling of references

    
    %%
    disp('Finding vesicles');
    
%     eMask=meGetMask(t.mi,size(ccmx));
    
    if displayOn
        %         subplot(2,3,2);
        %         imacs(ccmx);
        subplot(2,3,2);
        imacs(t.ccmx.*t.mask);
        drawnow;
    end;
    if findInMask
        t.ccmx2=t.ccmx;
    else
        t.ccmx2=t.ccmx.*t.mask;
    end;
    %     t.mask=eMask;
    
    if displayOn
        t.msub=t.ms;
    end;
    t.umodel=zeros(ns0);
    
    %%    % scan all the CC peaks
    t.nfound=0;
    t.mi.vesicle.x=[];
    t.mi.vesicle.y=[];
    t.mi.vesicle.r=[];
    t.mi.vesicle.s=[];
    t.mi.vesicle.ok=[];
    
else % m is a string
    switch lower(m)
        case 'next'
            maxN=mi;        % pick up the alternate arguments.
            thresh=rPars;
            ns=size(t.ccmx2,1);
            radScalings=t.rads.^-rExponent;  % undo scaling of references
            
            nf0=t.nfound;
            while t.nfound<nf0+maxN
                [t.globalmax, jx, jy]=max2d(t.ccmx2);
                if t.globalmax<(thresh(1)*t.rads(1)^rExponent)  % below minimum value
                    break
                end;
                
                %                 % find out which reference gave rise to the maximum
                %                 refi=t.ccmi(jx,jy);
                
                % Get the corrected peak value and interpolated radius
                % index
                [ampi, refii]=max1di(squeeze(t.ccs(jx,jy,:)).*radScalings);
% ampi
                refri=t.fitmin+(refii-1)*t.rstep;  % interpolated radius in pixels
                blank=fuzzymask(ns,2,blankRadiusFactor*refri+t.mbnThickness/2,t.mbnThickness,[jx jy]);
                t.ccmx2=t.ccmx2.*(1-blank);
% imacs(t.ccmx2);
% title(ampi);
% drawnow;
                if ampi<thresh(2) && ampi >thresh(1)*(1-maxFracMasked) && refii>1.5  % might be ok, check for overlap with mask
%                     Also we reject radii that are at the minimum.
                    refi=max(1,min(round(refii),t.nrsteps)); % closest model radius
                    testRef=circshift(t.orefs(:,:,refi),round([jx jy]-ns/2-1));
                    fracMasked=(testRef(:)'*(~t.mask(:)))/t.opowers(refi);
                    if fracMasked<maxFracMasked && ampi>thresh(1)*fracMasked
                        flag=single((fracMasked<.01) && (refri>=t.rmin)...
                         && (refri<=t.rmax) && ampi>thresh(1)/(1-fracMasked));  % in bounds: flag=1.
                        flag=flag & t.mask(jx,jy);
                        t.nfound=t.nfound+1;
                        t.mi.vesicle.r(t.nfound,1)=refri*t.ds;
                        t.mi.vesicle.x(t.nfound,1)=(jx-1)*t.ds+1;
                        t.mi.vesicle.y(t.nfound,1)=(jy-1)*t.ds+1;
                        t.mi.vesicle.s(t.nfound,1)=ampi/(1-fracMasked);
                        t.mi.vesicle.ok(t.nfound,1:4)=[1 flag flag flag];  % flag indicates a vesicle in bounds.
                        vref=ampi*circshift(t.rrefs(:,:,refi),round([jx jy]-ns/2-1));
                        t.umodel=t.umodel+Crop(vref,t.ns0);  % approximate model
                    end;
                end;
            end;
            if displayOn
                subplot(2,3,4); imacs(t.msub);
                title(t.nfound);
                subplot(2,3,2); imacs(t.umodel);
                subplot(2,3,6); imacs(t.ccmx2);
                
                subplot(2,3,3);
                plot(t.mi.vesicle.r*t.mi.pixA,t.mi.vesicle.s,'k.');
                xlabel('Vesicle radius, Å');
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
