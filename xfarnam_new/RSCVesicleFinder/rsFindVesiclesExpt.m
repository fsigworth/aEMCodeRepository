% function [mi, state]=rsFindVesiclesExpt(m,mi,rPars,findInMask)
% Auto vesicle finder for Vesicle_finding_GUI
% Calls:
% [mi t]=meFindVesicles3(m, mi, rPars, findInMask)          --initialize vesicle finder
% mi = meFindVesicles3('next',maxN,minmaxThresh);  --find up to maxN vesicles.
% mi = meFindVesicles3(m, mi) -- new image, old parameters
% meFindVesicles('end');      --deallocates the persistent variables.
% t is a persistent structure that contains information for the

findInMask=1;
displayOn=1;
rPars=[100 300 20];

basePath='/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/';
cd(basePath)
mi=ReadMiFile('Info/sq10_347_Jun25_22.44.49mi.txt');
mi.basePath=basePath;
m=meReadMergedImage(mi);
n0=[768 768];
m=Downsample(m,n0);

blankRadiusFactor=.7;
borderFactor=1.25;
maxFracMasked=.6;   % tries to fit vesicles up to this overlap with mask.
extendedFitRadius=1.6;
reducedFitRadius=0.8;
sHP=.01;  % highpass frequency in downsampled pix^-1
varLP=.02;  % lowpass for image var

targetPixA=10;
% displayOn=0;
rExponent=.2;  % radius-weighting

    t.mi=mi;
    m=GaussHP(m,sHP);
    
    % Get image and pixel sizes
    n0=size(m);
    ds0=mi.imageSize(1)/n0(1);  % downsampling factor of m
%     pixA0=mi.pixA*ds0;    % pixel size of m
    ns=n0; % use the original image size
%     ns0=NextNiceNumber(n0*pixA0/targetPixA);  % an even number, max about 10% too big.
    nsx=borderFactor*ns;  % expanded number of pixels
%     [t.msx, mag]=DownsampleGeneral(m,ns,ns0/n0);  % resample the image, downsample factor ds
    t.msx=Crop(m,nsx);  % pad the image
% t.ms=Crop(t.msx,ns0); % non-padded image
    t.ms=m;
    borderShift=(nsx-ns)/2;  % displacement of the center of padded image
    t.ns0=ns;
%     t.ds=ds0/mag;  % total downsampling (number > 1) from original image
    t.ds=ds0;
    pixA=t.ds*mi.pixA;  %pixA in the image ms
    
    vm=meDownsampleVesicleModel(mi.vesicleModel,t.ds)*pixA;
    
    %   note that the density is scaled up to match the scattering per voxel.
    t.mbnThickness=(sum(vm>max(vm)/2));  % approx. membrane thickness in pixels.
    
%     Make a Gaussian HP filter
    highPass=exp(-sHP^2./(max(RadiusNorm(nsx),1e-6).^2));  % zero-center highpass function
    t.fms=fftn(t.msx).*ifftshift(highPass);  % Fourier-transformed, filtered image.
    
    if displayOn
        SetGrayscale;
        subplot(2,3,1);
        imacs(t.ms);
    end;
    
    % Get the effective CTF from the merging.
    H=ifftshift(highPass.*meGetEffectiveCTF(mi,nsx,t.ds));
    
    % Get the mask, pad it to the search size
    t.mask=zeros(nsx);
    t.mask(1:ns(1),1:ns(2))=meGetMask(mi,ns);
    t.findInMask=findInMask;
    
    %%
    
%     if ~(isfield(t,'rPars') && all(t.rPars==rPars) ...
%             && nsx(1) == size(t.frefs,1))  % need to compute references
        disp('Constructing references...');
        t.rPars=rPars;
        t.rmin=rPars(1)/pixA;  % minimum liposome radius in pixels
        t.rmax=rPars(2)/pixA;  % maximum liposome radius
        t.fitmin=t.rmin*reducedFitRadius;
        t.fitmax=t.rmax*extendedFitRadius;
        t.rstep=rPars(3)/pixA; % radius step
        t.nrsteps=round((t.fitmax-t.fitmin)/t.rstep+1);
        
        t.frefs=single(zeros([nsx,t.nrsteps]));
        t.localSD=zeros([ns,t.nrsteps],'single');
        t.rads=zeros(t.nrsteps,1);
        t.powers=zeros(t.nrsteps,1);
        t.rrefs=single(zeros([nsx t.nrsteps]));
        t.orefs=single(zeros([nsx t.nrsteps]));
        t.opowers=zeros(t.nrsteps,1);
        
        if displayOn
            
            subplot(2,3,2);
            nd=NextNiceNumber(3*t.rmax);  % size of box to display
            nd=min(nd,nsx);
        end;
        vesOrigin=(nsx/2+1)+borderShift;   % shift the origin of the references
%         so they'll correspond to the center of the padded image.  zero
%         lag will then be (1,1), pointing to the correct location of the
%         unpadded image.
        mVar=GaussFilt(t.ms,varLP).^2;
        fmVar=fftn(mVar);
        for i=1:t.nrsteps
%             Make vesicle model
            r=t.fitmin+(i-1)*t.rstep;  % radius in pixels
            t.rads(i)=r;  % radius in A
            v = VesicleFromModel(nsx,r,vm,vesOrigin);
            
%             Compute the local SD, ns points
            rd=ifftshift(fuzzymask(ns,2,r+1.5*t.mbnThickness,t.mbnThickness/2));
            fd=fftn(rd); % disc and ft for variance calc
            lsd=sqrt(abs(real(ifftn(fmVar.*fd))/(rd(:)'*rd(:))));
            t.localSD(:,:,i)=lsd;

            fv=-H.*fftn(ifftshift(v));    % FT of vesicle at origin.
            t.frefs(:,:,i)=fv;
            t.orefs(:,:,i)=v;
            t.opowers(i)=v(:)'*v(:);
            rv=fftshift(real(ifftn(fv)));  % CTF-filtered real image of vesicle, centered
            t.powers(i)=rv(:)'*rv(:);
            t.rrefs(:,:,i)=rv;
            
            if displayOn && mod(i,10)==0
                %                 t.rrefs(:,:,i)=rv;
                imacs(rv);
                title(i);
                drawnow;
            end;
        end;
%     end;
    medSD=median(t.localSD(:));
%     Make localSD bigger.
    t.localSD(ns(1)+1:nsx(1),:,:)=1;
    t.localSD(:,ns(2)+1:nsx(2),:)=1;
    %%
    % Compute all the cross-correlations, using a weight of r^rExponent to equalize
    % the CC peak values
    disp('Computing cross-correlations');
    t.ccs=single(zeros([nsx t.nrsteps]));
    t.nccs=single(zeros([nsx t.nrsteps]));
    % cc values are scaled by 1/r^rExponent to make small vesicles as easily found
    for i=1:t.nrsteps
        t.ccs(:,:,i)=t.rads(i).^(rExponent)/t.powers(i)*...
            real(ifftn(t.fms.*conj(t.frefs(:,:,i))));
        t.nccs(:,:,i)=t.ccs(:,:,i)./t.localSD(:,:,i)*t.rads(i)^-rExponent;
    end;
    % ccs=ccs-min(ccs(:));  % min value is zero.
    [t.ccmx, t.ccmi]=max(t.ccs,[],3);
    t.ccmxScaled=t.ccmx.*(t.rads(t.ccmi).^-rExponent);  % undo scaling of references

    [t.nccmx, t.nccmi]=max(t.nccs,[],3);
t.nccmx=t.nccmx*medSD;
% t.ccmx=t.nccmx;    
    %%
    disp('Finding vesicles');
    
%     eMask=meGetMask(t.mi,size(ccmx));
    
    if displayOn
        %         subplot(2,3,2);
        %         imacs(ccmx);
        subplot(2,3,2);
        imacs(t.ccmx.*t.mask);
        subplot(2,3,3);
        imacs(t.nccmx);
        subplot(2,3,4);
        imacs(mVar);
        subplot(2,3,5);

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
    t.umodel=zeros(ns);
    
    %%    % scan all the CC peaks
    t.nfound=0;
    t.mi.vesicle.x=[];
    t.mi.vesicle.y=[];
    t.mi.vesicle.r=[];
    t.mi.vesicle.s=[];
    t.mi.vesicle.ok=[];

    
    return
% else % m is a string, 'next' or 'end'
    switch lower(m)
        case 'next'
            maxN=mi;        % pick up the alternate arguments.
            thresh=rPars;
            nsx=size(t.ccmx2,1);
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
ampi=t.globalmax; %%%%%%
% ampi
                refri=t.fitmin+(refii-1)*t.rstep;  % interpolated radius in pixels
                blank=fuzzymask(nsx,2,blankRadiusFactor*refri+t.mbnThickness/2,t.mbnThickness,[jx jy]);
                t.ccmx2=t.ccmx2.*(1-blank);
% imacs(t.ccmx2);
% title(ampi);
% drawnow;
                if ampi<thresh(2) && ampi >thresh(1)*(1-maxFracMasked) && refii>1.5  % might be ok, check for overlap with mask
%                     Also we reject radii that are at the minimum.
                    refi=max(1,min(round(refii),t.nrsteps)); % closest model radius
                    testRef=circshift(t.orefs(:,:,refi),round([jx jy]-nsx/2-1));
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
                        vref=ampi*circshift(t.rrefs(:,:,refi),round([jx jy]-nsx/2-1));
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
                plot(t.mi.vesicle.r(:,1)*t.mi.pixA,t.mi.vesicle.s(:,1),'k.');
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
