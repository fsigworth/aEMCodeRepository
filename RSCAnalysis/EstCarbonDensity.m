% function [s,sp1,=EstCarbonDensity(mi)
% function s=EstCarbonDensity(mi)
% Given the mi structure, read a merged image, subtract vesicles, compute
% 1D spectrum, fit the carbon spectrum+shot noise

% [stackName,stackPath]=uigetfile('*si.mat');
% [basePath,stackDir]=ParsePath(stackPath);
% cd(basePath);
% 
% load([stackDir stackName]);

useNonParticles=0;

fc1=.05;
    fc2=.01;
    imgFraction=.7;
s=struct;

if useNonParticles
[si,imgs,fnames,siPath]=reLoadStackFiles;
n1=size(imgs,1);
nmi=numel(si.mi);
disp([num2str(nmi) ' mi files.']);
end;
%%
startInd=2;
lastInd=2;
for ind=startInd:lastInd
    ind
    mi=si.mi{ind};
if useNonParticles    
%     find the stack images corresponding to this image, and compute their
%     spectrum.  Note that particle images are already scaled up by
%     sqrt(doses(1))*pixA/ds
    partFlags=(si.miIndex==ind);
    
    sp1n=RadialPowerSpectrum(imgs(:,:,partFlags),1);
    sp1m=mean(sp1n,2);
    freqs1=(0:n1/2-1)/(n1*si.pixA);
    subplot(2,2,1);
    plot(freqs1,sp1n);
else
    n1=128;
    freqs1=(0:n1/2-1)/(n1*si.pixA);
    sp1n=0*freqs1';
    sp1m=0*freqs1';
    subplot(2,2,1);
    cla;
end;
    %%
    disp('Reading the merged image.');
    m=meReadMergedImage(mi);
    n=size(m);
    ds=mi.imageSize(1)/n(1);
    subplot(2,2,2);
    imags(GaussFilt(m,.1));
    drawnow;
    
    disp('Making vesicles');
    mv=m-meMakeModelVesicles(mi,m);
    imags(mv);
    drawnow;
    %%
    disp('Making mask');
    var=GaussFilt(GaussFilt(mv,fc1).^2,fc2);
    thresh=Percentile(var,imgFraction);
    msk=GaussFilt(GaussFilt(var<thresh,fc2)>.95,fc2);
    subplot(2,2,1+useNonParticles);
    imags(mv.*msk);
    s0=median(mi.vesicle.s(:,1));
        title([num2str(ind) ':   s0= ' num2str(s0,3) '   dose=' num2str(mi.doses(1),3)]);

    drawnow;
    
    sp2=fftshift(abs(fft2(msk.*mv)).^2)*mi.doses(1)/(msk(:)'*msk(:)*mi.pixA*ds);
    [sp1,freqs]=RadialCTF(sp2,mi.ctf(1),mi.pixA*ds);
    
    mi1=mi;
        mi1.ctf(1).B=30;
    % Remove astigmatism for 1D model
    for i=1:numel(mi.ctf)
        mi1.ctf(i).deltadef=0;
    end;
    %
    mi1.weights(2)=1;
    ec=(meGetEffectiveCTF(mi1,n)).^2;
    %
    
    n2=numel(freqs);
    nf0=find(freqs>.03,1);  % min frequency for fitting
    nf1=round(n(1)*.8/2);  % fit up to 80% of nyquist
    nf=nf1-nf0+1;
    pars=[1 .031 3.3];
    pars=[1 .031 6];
    
    F=ones(n2,4);
    F(:,1)=(0:n2-1)'/n2;
    F(:,3)=sectr(ec);
    cs=CarbonSpectrum(freqs,pars);
    F(:,4)=sectr(ec).*(cs-1);
    a=LinLeastSquares(F(nf0:nf1,:),sp1(nf0:nf1));
    
    model=F*a;
    subplot(2,1,2);
    mxv=max(model);
    plot(freqs,min(1.2*mxv,[sp1 model a(3)+a(4)*(cs-1)]));
    hold on;
    if useNonParticles
        plot(freqs1,mean(sp1n,2));
    end;
    plot([freqs(nf0) freqs(nf1)],[.5 .5]);
    hold off;
    title([num2str(ind) ' fit:  ' num2str(a(2:4)',3)...
        '   s0= ' num2str(s0,3)...
     '   dose=' num2str(mi.doses(1),3) '  exp= ' num2str(pars(3))]);
    drawnow;
    
        subplot(2,2,2);
    imags(GaussFilt(m,.1));

    
    s.a=a';
    s.carbon=pars;
    mi.carbonModel=s;
    si.mi{ind}=mi;
    disp('done.');
end;
return
