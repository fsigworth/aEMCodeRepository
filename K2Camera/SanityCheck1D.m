% SanityCheck1D
nIter=1000;
err=zeros(nIter,1);
dx=zeros(nIter,1);

for iter=1:nIter
    
    n=1024;
    sigAmp=1;
    hw=30;
    noiseAmp=.05;
    fcn=.05;
    ctr=n/2+1;
    t=sect(Radius(n));
    f=t/n;
    % motif=max((t-n/2-1+hw),0);
    % derivative is sigAmp * 1 unit/pixel
    
    sigma=10/n;
    delta=1/sqrt(2*pi)*exp(-f.^2/(2*sigma^2));
    motifCtr=GaussFilt(delta,.02);
    motif=circshift(motifCtr,n/2);
    
    % deriv1=diff(motifCtr);
    % deriv2=diff(deriv1);
    % val=motifCtr(ctr);
    % d1=deriv1(ctr)
    % d2=deriv2(ctr)
    
    n1=GaussFilt(randn(n,1),fcn)*noiseAmp;
    n2=GaussFilt(randn(n,1),fcn)*noiseAmp;
    i1=n1+sigAmp*motif;
    i2=n2+sigAmp*motif;
    
    cs = fft(i1).*conj(fft(i2));
    css=fftshift(cs);
    r=real(fftshift(ifft(cs)));  % ifft has an implicit divide by n
    subplot(231);
    plot(Crop(r,17));
    [val x]=max1di(r);
    dx(iter)=x-ctr;
    
%     disp('true derivatives');
    deriv1=diff(r);
    deriv2=diff(deriv1);
    val=r(ctr);
    d1=deriv1(ctr);
    d2=deriv2(ctr);
    
    
    
    
    
    f=t/n;
    
    % subplot(231); plot(r);
    subplot(232); plot([real(css) imag(css)]);
    %%
    %
%     disp('computed');
    d1=2*pi*sum(real(css).*f)/n;
    d2=-4*pi^2*sum(real(css).*f.^2)/n;
    sd1=4*pi^2*sum(2*imag(css).^2.*f.^2)/n^2;
    err(iter)=(sd1/d2^2);
    
end;

std(dx)
mean(sqrt(err))

