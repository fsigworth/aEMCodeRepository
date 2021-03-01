function img=rspReconstructParticles(dis,mi,picks,ptrs,rscc)
% For SimpleRSPicker, create an image of autopicked particles, based on the
% best-matching template and the particle amplitudes.
scaleUp=2;
minH=.25;  % specify the maximum magnification at low frequencies
if isfield(rscc,'fHP')
    fHP=rscc.fHP;
else
    fHP=.001; % highpass in A^-1 used in picking preprocessor
end;

img=single(zeros(dis.ndis));
if ~isfield(rscc,'eigenImgs') || ~isfield(rscc,'vList') % data not there
    return
end;
apFlags=16:47;  % Flag for picked particles
apPtrs=[2 3 4 6];    % pointers for picked particles
n=size(rscc.mxCC);
ds1=mi.imageSize(1)/n(1);
[ne, ~, nterms]=size(rscc.eigenImgs);
vsz=numel(rscc.vList)/nterms;
eigs=reshape(rscc.eigenImgs,ne*ne,nterms);
vlst=reshape(rscc.vList,nterms,vsz);

scaled=0;
% Create all the model particles before filter compensation.
for ptrIndex=apPtrs
    for ip=1:ptrs(ptrIndex)  % loop over all auto-picked particles
        c=squeeze(picks(ptrIndex,ip,:))';
        pos=c(1:2)/ds1+1;
        flag=c(3);
        vesInd=c(4);
        s=0;
        if vesInd>0
            s=mi.vesicle.s(vesInd,1);
        elseif dis.spMode
            s=1;
        end;
        if dis.useRawAmplitudes
            amp=-c(5);              % invert the contrast
            ampU=amp/dis.ccuScale;  % unscaled (absolute) particle amplitude
        else
            amp=-c(5);
            ampU=amp*s;
        end;
        templ=c(6);  % template index
        if any(flag==apFlags) && templ>0 && templ<=size(vlst,2)
            
            % We could call rspMakeModelParticle here, but probably not as
            % efficient.
            modelParticle=reshape(eigs*vlst(:,templ),ne,ne);
            %             if numel(dis.invPWF)>1  % undo the prewhitening filter
            %                 modelParticle=real(ifftn(fftn(modelParticle).*ifftshift(dis.invPWF)));
            %             end;
            img=img+ampU*scaleUp*ExtractImage(modelParticle,round(pos),n,1)/2;
        end;
    end;
end;
if scaled
    disp('Scaled amplitudes used');
end;

% Undo the prewhitening (and highpass) filter over the entire image.
if isfield(rscc,'pwImg')
    h=rscc.pwImg;
else
    h=GetNoiseWhitening(mi,dis.ndis,dis.ds,0);
    h(h<minH)=minH; % we don't want to boost too much in the inverse filter.
end;
img=real(ifftn(fftn(img)./ifftshift(h)));

% do the inverse-CTF boost, and display lowpass and highpass
ctComp=dis.filter(3)/100;  % inverse CTF amplitude, given in percent
fc=mi.pixA*dis.ds./dis.filter; % Convert from A to Fourier pixel
fc(dis.filter==0)=0;           % 0 A -> 0 A^-1

% The CTFInverseFilter also does the highpass at fc(1).
img=rspCTFInverseFilter(img,mi,ctComp);
% fLP=fc(2)*dis.ds*mi.pixA; % don't lowpass the model
% img=GaussFilt(img,fc(2)); % lowpass like the display
img=img*dis.mulr+dis.addr;  % scale the samee as the image.
% dis.mode=3;
% imgs(:,:,3)=img;
%             rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);

% Copied from rsPickingPreprocessor4
function T=GetNoiseWhitening(mi,np,ds,fHP)
pwFiltPars=[.002 0; .01 .5];  % generic pw filter parameters
T=meGetNoiseWhiteningFilter(mi,np,ds,1,fHP*mi.pixA*ds);
if numel(T)<2
    disp('Using a generic PW filter function.');
    freqs=RadiusNorm(np)/mi.pixA+1e-6;  % in inverse Å, prevent divide by zero
    T=ones(np,'single');
    for i=1:size(pwFiltPars,1)  % product of (gauss + const) terms.
        f0=pwFiltPars(i,1);
        a=pwFiltPars(i,2);
        h=exp(-(f0./freqs).^2);
        T=T.*(a+(1-a)*h);
    end;
end;

