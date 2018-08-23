% SanityCheck3
nIter=1;
nj=10;
n=1024;
ndis=64;
dr0=zeros(nIter,1);
drw=zeros(nIter,1);
err=zeros(nIter,1);
fc=80;
for iter=1:nIter
    sigAmp=1;
    hw=30;
    noiseAmp=.4;
    fcn=.1;
    ctr=n/2+1;
    t=Radius(n);
    f=t/n;
    % motif=max((t-n/2-1+hw),0);
    % derivative is sigAmp * 1 unit/pixel
    
    sigma=10/n;
    delta=1/(2*pi)*exp(-f.^2/(2*sigma^2));
    motifCtr=GaussFilt(delta,.02);
    motif=ifftshift(motifCtr);
    imgs=zeros(n,n,nj);
    imgMeans=zeros(n,n,nj);
    for j=1:nj
        imgs(:,:,j)=GaussFilt(randn(n),fcn)*noiseAmp;
    end;
    n2=GaussFilt(randn(n),fcn)*noiseAmp;
    imgs=imgs*noiseAmp+sigAmp*repmat(motif,[1 1 nj]);
    
    sumImg=sum(imgs,3);

    
    
    meanR=zeros(n,n);
    meanI=zeros(n,n);
    meanC=zeros(n,n);
    for j=1:nj
        meanImg=(sumImg-imgs(:,:,j))/(nj-1);
        imgMeans(:,:,j)=meanImg;
        cs=fftn(imgs(:,:,j)).*conj(fftn(meanImg));
        meanR=meanR+real(cs);
        meanI=meanI+imag(cs).^2;
        meanC=meanC+cs;
    end;
    
[H1, rs]=k2GetWeightingFilter3(imgs,imgMeans,0);
    return
    
    
    meanC=meanC/nj;
    meanR=meanR/(nj);
    meanR1=Radial(fftshift(meanR));
    meanI=meanI/(nj);
    meanI1=Radial(fftshift(meanI));
    
    w=f1.*(meanR1.*f1*2*pi/n)./(meanI1.*f1*2*pi/n^3);
    w(1)=0;  % otherwise NaN
    subplot(234)
    % w1=Radial(w);
    w(fc:512)=0;
    w=GaussFilt(w,.1);
    w=max(w,0);
    plot(w(1:ndis));
    title('Filter W');
    
    %
    %     cs = fftn(i1).*conj(fftn(i2));
    %     css=fftshift(cs);
    r=real(fftshift(ifftn(meanC)));  % ifft has an implicit divide by n
    subplot(231);
    imacs(Crop(r,n));
    [val, x, y]=max2di(r);
    dx=[x y]-ctr;
    rawShift=hypot(x-ctr,y-ctr)
    dr0(iter)=rawShift;
    
    
    
    
    d2rw=8*pi^3*cumsum(meanR1.*f1.^3.*w)/n;
    
    var1rw=8*pi^3*sum(meanI1.*f1.^3.*w)/n^3
    
    var2rw=16*pi^4*cumsum(meanI1.*f1.^5.*w.^2./n^2);
    
    sigma2=d2rw./sqrt(var2rw);
    subplot(236);
    plot(sigma2(1:ndis));
    title('sigma2');
    
    sdW=sqrt(var1rw)/d2rw;
    
    Lnum=8*pi^3*cumsum((meanR1.*f1.^3.*w).^2);
    Lden=8*pi^3*cumsum(meanI1.*f1.^3.*w)/n^3;
    L=sqrt(Lnum./(Lden*sqrt(n)) );
    subplot(235)
    plot(L(1:ndis));
    title('Precision l');
    w2=ToRect(w);
    Lmax=L(ndis)
    rw=real(fftshift(ifftn(cs.*ifftshift(w2))));  % ifft has an implicit divide by n
    subplot(232);
    imacs(Crop(rw,n));
    rw=rw*(max(r)/max(rw));
    subplot(233)
    plot([sect(Crop(rw,128)) sect(Crop(r,128))]);
    [val, x, y]=max2di(rw);
    dxw=[x y]-ctr;
    rawShiftW=hypot(x-ctr,y-ctr)
    drw(iter)=rawShift;
    
end;   
    
    return
    %
    %     %     disp('true derivatives');
    %     [d1y,d1x]=gradient(r);
    %     [dxy, dxx]=gradient(d1x);
    %     val=r(ctr,ctr);
    %     d1=d1x(ctr,ctr);
    %     d2=dxx(ctr,ctr);
    %
    %     %     subplot(231);
    %     %     imacs(r);
    %     subplot(232)
    %     plot(sect(Crop(d1x,128)));
    %     subplot(233);
    %     plot(sect(Crop(dxx,128)));
    %
    %
    %     f=t/n;
    %     fx=repmat(sect(f),1,n);
    %     % subplot(231); plot(r);
    %     subplot(232); plot([real(css) imag(css)]);
    %%
    %
    %     disp('computed');
    %     d1x=2*pi*sum(real(css).*f)/n^2;
    %     d2=4*pi^2*sum(sum(real(css).*f.^2))/n^2
    %     var1=4*pi^2*sum(sum(2*imag(css).^2.*f.^2))/n^4
    %     sd1=sqrt(var1);
    %     err(iter)=var1/d2^2;
    %
    %     f1=sectr(f);
    %     %     circ sums
    %     d2r=8*pi^3*sum(Radial(real(css)).*f1.^3)/n
    %
    %     var1r=8*pi^3*sum(Radial(2*imag(css).^2).*f1.^3)/n^3
    %
    %     w=f1.*(Radial(real(css)).*f1*2*pi/n)./sqrt(Radial(2*imag(css).^2).*f1*2*pi/n^3);
    %     w(1)=0;  % otherwise NaN
    %     subplot(234)
    %     % w1=Radial(w);
    %     w(fc:512)=0;
    %     w=GaussFilt(w,.1);
    %     w=max(w,0);
    %     plot(w(1:ndis));
    %     title('Filter W');
    %
    %     d2rw=8*pi^3*sum(Radial(real(css)).*f1.^3.*w)/n
    %
    %     var1rw=8*pi^3*sum(Radial(2*imag(css).^2).*f1.^3.*w)/n^3
    %
    %     sdW=sqrt(var1rw)/d2rw
    %
    %     Lnum=8*pi^3*cumsum((Radial(real(css)).*f1.^3.*w).^2);
    %     Lden=8*pi^3*cumsum(Radial(2*imag(css).^2).*f1.^3.*w)/n^3;
    %     L=sqrt(Lnum./(Lden*sqrt(n)) );
    %     subplot(235)
    %     plot(L(1:ndis));
    %     title('Precision L');
    %     w2=ToRect(w);
    %
    %     rw=real(fftshift(ifftn(cs.*w2)));  % ifft has an implicit divide by n
    %
    %     subplot(236)
    %     plot(sect(Crop(rw,128)));
    %
    %
    %
    % end;
    % measSD=std(dr)
    predSD=sqrt(mean(err))  % about sqrt(2) too big.