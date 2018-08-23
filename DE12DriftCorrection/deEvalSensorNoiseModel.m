function specModel=deEvalSensorNoiseModel(s,ds)
% function specModel=deEvalSensorNoiseModel(s,ds)
%  obtain the **zero-centered** sensor noise spectral model for images downsampled by the
%  factor ds.  That is, a corrected cross-spectrum is obtained from images
%  i1 and i2 of size nw and scaled the same as individual frames in s,
%    corrCrossSpectrum=ifftn(fftn(i1).*conj(fftn(i2)))/prod(nw)-ifftshift(specModel);

    [nx ny nim]=size(s.imgs);
    nxd=nx/ds;  % downsampled image size
    nyd=ny/ds;
    p=s.modelPars;
    np=numel(s.modelPars);
    
    bfmat=zeros(nxd*nyd,np);
 
    [x y]=ndgrid((-nxd/2:nxd/2-1)'/nx,(-nyd/2:nyd/2-1)'/ny);
    
    for j=1:np
        if j<np fvals=x; else fvals=y; end;
        f=(exp(-abs(p(j)*fvals).^2/2));
        bfmat(:,j)=f(:);
    end;
%     Correct the cropped cross-spectrum to match the effect of downsampling
    specModel=reshape(bfmat*s.modelAmps,nxd,nyd)/ds^2;
