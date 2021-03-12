% HoleEdgeMasker.m

doWrite=1;


allMisName='Picking_9/allMis9_intens+frac_7505.mat';
outMisName='Picking_9/allMis_holemasked_i+f_7505.mat';
outJPath='Picking_9/MaskJpegs/';
outLocsName='Picking_9/MaskLocs.mat';

% % disp(['Loading ' allMisName]);
% % load(allMisName);
nmi=numel(allMis)
outMis=cell(nmi,1);
%%
fHP=0;
k=.1;
priorAmp=100;
if doWrite
       CheckAndMakeDir(outJPath,1);
end;

dds=8;
r=1;
nBig=768;

% priorPos=[350 380 219 264 507 545
%     370 170 337 555 554 344];
% priorSigma=50;

%
nHoles=1;
locs=zeros(ceil(nmi/6),2,6,'single');
%
update2DFreq=5000

figure(1);
iStart=1
for i=iStart:numel(allMis)
    mi=allMis{i};
    imgName=['Merged_sms/' mi.baseFilename 'ms.mrc'];
    disp([imgName '  ' num2str([nHoles i])]);
    if ~exist(imgName,'file')
        disp(' no image.');
        continue;
    end;
    iShot=str2double (mi.baseFilename(end))+1;
    
    % load(['/Volumes/D257/Hideki/201228/05_1_1/Merged_sm/' baseFilename 'rscc.mat']);
    m0=ReadMRC(imgName);
    % mi=ReadMiFile(['Info_4_pikPart1/' baseFilename 'mi.txt']);
    imags(m0);
    
    m0=m0-mean(m0(:));
    %     m0=GaussHP(m0,fHP);
    n0=size(m0);
    ds=mi.padImageSize(1)/n0(1);
    m0PixA=mi.pixA*ds;
    n8=n0/8;
    m8=Downsample(m0,n8);
    
    if i==1
        pixA8=8*m0PixA;
        
        rPix=r*1e4/pixA8;
        ref=fuzzymask(nBig,2,rPix,0);
        nccMask=fuzzymask(nBig,2,rPix,1);
        % ref=GaussHP(ref,fHP);
        box=ones(n8,'single');
        boxp=Crop(box,nBig);
        fref=fftn(ref);
        nc=real(ifftn(conj(fftn(ifftshift(boxp))).*fref));
        % subplot(232);
        % imags(cc);
        
        p=nc/prod(n8);
        norm=p.*(1-p);
        
        subplot(231);
        imags(ref);
        subplot(232);
        imags(1./(norm+k*10));
        plot((1./(norm(:,1:20:nBig)+k)));
        
    end;
    
    
    % create the filters for detecting carbon.
    fHP=.2; % in A^-1
    kHP=-(log(2)/2)*(m0PixA*fHP)^2;
    R=RadiusNorm(n0);
    q=exp(kHP./(R+1e-4));	% Gaussian kernel
    cPars=mi.ctf;
    cPars.B=40;
    cPars.alpha=.05;
    c1=CTF(n0,m0PixA,cPars).*q;
    cPars.alpha=pi/2+.05;
    c2=CTF(n0,m0PixA,cPars).*q;
    
    fm0=fftshift(fftn(m0)); % zero center FT
    mf1=real(ifftn(ifftshift(c1.*fm0)));
    mf2=real(ifftn(ifftshift(c2.*fm0)));
    
    % %     figure(3);
    %      sp1=RadialPowerSpectrum(fm1);
    % sp2=RadialPowerSpectrum(fm2);
    % freqs=sectr(RadiusNorm(n0))/m0PixA;
    % semilogy(freqs,[sp1 sp2]);
    
    %  Compute the local variance
    md1=max(0,Downsample(GaussFilt(mf1.^2,.05),n8));
    md2=max(0,Downsample(GaussFilt(mf2.^2,.05),n8));
    
    v8=sqrt(md2)-sqrt(md1);
    v8=v8-median(v8(:));
    % plot([sectr(16*v8) sectr(m8)]);
    % Combine variance and mean image darkness
    varWeight=2*16; % about twice contribution from var
    mv8=m8+(varWeight*v8);
    
    m8p=Crop(mv8,nBig);
    % m8p=GaussHP(m8p,fHP);
    
    subplot(234);
    imags(m8);
    subplot(235);
    imags(mv8);
    cc=real(ifftn(conj(fftn(ifftshift(m8p))).*(fref)));
    % subplot(232);
    % imags(cc);
    ncc0=cc./(norm+k);
    
    %    prior=priorAmp*Gaussian(nBig,2,priorSigma,priorPos(:,iShot)');
    
%     subplot(232);
%     plot([sect(ncc0) sect(ncc0')]);
%     n81=n8(1);
%     axis([1 nBig -n81 n81]);
    
    prior=priorAmp*Gaussian(nBig,2,nBig/4);
    ncc=(ncc0+prior).*nccMask;
     
    [mx,ix,iy]=max2d(ncc);
    
    subplot(233);
    imags(ncc);
    hold on;
    plot(ix,iy,'r+','markersize',10);
    hold off;
    
    model=ExtractImage(ref,[ix,iy],n8);
    subplot(236);
    
    modelx=GaussFiltDCT(model,.05)>.9; % expand by a few pixels
    outMis{i}=meInsertMask(modelx,mi,2);
    
    if doWrite
        jname=sprintf('%s%05d.jpg',outJPath,i);
        modelImage=Downsample(modelx+60*m8,4*n8);
        WriteJpeg(modelImage,jname);
        disp(jname);
        
        %         mi=meInsertMask(modelx,mi,2);
        %         WriteMiFile(mi,[outMiPath mi.baseFilename 'mi.txt']);
    end;
    
    imags(modelx+30*mv8);
    
    locs(nHoles,:,iShot)=[ix iy];
    
    title(iShot);
    % pause;
    
    if iShot==1
        nHoles=nHoles+1;
    end;
    if mod(nHoles,update2DFreq)==0
        xs=squeeze(locs(1:nHoles,1,:));
        ys=squeeze(locs(1:nHoles,2,:));
        figure(2);
        plot(xs,ys,'.');
        figure(1);
    end;
    drawnow;
end;
if doWrite
    disp(['Writing ' outLocsName '...']);
    save(outLocsName,'locs');
    disp(['Writing ' outMisName '...']);
    oldAllMis=allMis;
    allMis=outMis;
    save(outMisName,'allMis');
    allMis=oldAllMis;
    clear oldAllMis
    disp('done.');
end;
return


















%
% % save locs001amp60.mat locs
% locs=locs(1:nHoles,:,:);
% load locs001amp60.mat

xs=squeeze(locs(:,1,:));
ys=squeeze(locs(:,2,:));



figure(2);
plot(xs,ys,'.');
medx0=median(xs);
% medx=mean(squeeze(locs(:,1,:)));
medy0=median(ys);
% medy=mean(squeeze(locs(:,2,:)));
hold on;
plot(medx0,medy0,'ko');
hold off;
legend({'0' '1' '2' '3' '4' '5'},'fontsize',12);

figure(3);
shots=2:6;
xs1=xs(:,shots);
ys1=ys(:,shots);
subplot(211);
plot(xs1,'o');
medx=median(xs1);
hold on;
plot([1 nHoles],[medx;medx]);
hold off;
% legend({'0' '1' '2' '3' '4' '5'},'fontsize',12);
legend({'1' '2' '3' '4' '5'},'fontsize',12);
title('X coords');

subplot(212);
plot(ys1,'o');
medy=median(ys1);
hold on;
plot([1 nHoles],[medy;medy]);
hold off;
% legend({'0' '1' '2' '3' '4' '5'},'fontsize',12);
legend({'1' '2' '3' '4' '5'},'fontsize',12);
title('Y coords');

figure(2);
hold on;
plot(medx,medy,'ko');
hold off;

shots=1:6;
figure(4);
subplot(211);
hist(xs,100);
hold on;
plot(medx,500,'ko');
hold off;
xlabel('X values');
legend({'0' '1' '2' '3' '4' '5'},'fontsize',12);

subplot(212);
hist(ys,100);
hold on;
plot(medy,500,'ko');
hold off;
xlabel('Y values');

pause(1);

% Try editing values outside the inner ring of the reference
ctr=nBig/2+1;
bad=sqrt((xs-ctr).^2+(ys-ctr).^2)>rPix;
xs(bad)=NaN;
ys(bad)=NaN;
for i=1:6
    y=ys(:,i);
    x=xs(:,i);
    medy(i)=median(y(~isnan(y)));
    medx(i)=median(x(~isnan(x)));
end;

figure(2);
plot(xs,ys,'.');
hold on;
plot(medx,medy,'ko');
hold off;
legend({'0' '1' '2' '3' '4' '5'},'fontsize',12);

figure(3);
shots=2:6;
xs1=xs(:,shots);
ys1=ys(:,shots);
subplot(211);
plot(xs1,'o');
medx1=medx(shots);
hold on;
plot([1 nHoles],[medx1;medx1]);
hold off;
% legend({'0' '1' '2' '3' '4' '5'},'fontsize',12);
legend({'1' '2' '3' '4' '5'},'fontsize',12);
title('X coords');

subplot(212);
plot(ys1,'o');
medy1=medy(shots);
hold on;
plot([1 nHoles],[medy1;medy1]);
hold off;
% legend({'0' '1' '2' '3' '4' '5'},'fontsize',12);
legend({'1' '2' '3' '4' '5'},'fontsize',12);
title('Y coords');


figure(4);
subplot(211);
hist(xs,100);
hold on;
plot(medx,500,'ko');
hold off;
xlabel('X values');
legend({'0' '1' '2' '3' '4' '5'},'fontsize',12);

subplot(212);
hist(ys,100);
hold on;
plot(medy,500,'ko');
hold off;
xlabel('Y values');
