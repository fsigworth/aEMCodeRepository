function filter2D=hrpMakePWFilter(m1,pixA,displayOn)
% Create the pw filter for a micrograph
% -- currently works only for square images, need to fix ToRect.

if nargin<3
    displayOn=0;
end

n1=size(m1);

% Find local var
fc=.01*pixA;  % .01 A^-1 corner frequency for variance calc;
m1hp=GaussHP(m1,fc/10);
m1bp=GaussFilt(m1hp,fc);
localVar=GaussFilt(m1bp.^2,fc/4);
sz=size(m1bp);
msk=localVar<median(localVar(:));
msk1=GaussFilt(GaussFilt(msk,fc/2)>.9,fc/2);
% mysubplot(224);
% imags(msk1.*sqrt(localVar));
unmaskedFraction=sum(msk1(:).^2)/numel(m1);

if displayOn
    clf;
    mysubplot(221);
    imags(GaussFilt(m1,.1))
    mysubplot(222);
    imags(sqrt(localVar));
    mysubplot(223)
    imags(GaussFilt(m1hp.*msk1,.1));
    title(['Unmasked fraction = ' num2str(unmaskedFraction,3)]);
    drawnow;
end;

sp0=RadialPowerSpectrum(m1hp);
sp1=RadialPowerSpectrum(msk1.*m1hp)/unmaskedFraction;
ns=numel(sp0);
freqs=(1:ns)/(2*ns*pixA);

% semilogy(freqs,[sp0 sp1],'.-');

hfLimit=floor(2*ns*0.44); % rolloff from motionCor2 I guess.

pwFilter1=zeros(ns,1,'single');
pwFilter1(1:hfLimit)=1./sqrt(sp1(1:hfLimit));
pwFilter1(hfLimit+1:end)=pwFilter1(hfLimit);
fcDct=100/ns; % Fourier smoothing
pwFilter1=GaussFiltDCT(pwFilter1,fcDct);
filter2D=Crop(ToRect(pwFilter1),n1);

if displayOn
    mysubplot(224);
    plot(freqs,pwFilter1)

    mysubplot(222);
    imags(filter2D);
end;
end

