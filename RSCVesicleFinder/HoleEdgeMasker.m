% HoleEdgeMasker.m

fHP=0;
k=.001;
priorAmp=40;
runStr='_hole_all_amp20py4msk';
outMiPath=['Info' runStr '/'];
CheckAndMakeDir(outMiPath,1);
for i=1:6
    outJPaths{i}=sprintf('%s%s%s%1d/','jpeg',runStr,'shot',i);
    CheckAndMakeDir(outJPaths{i},1);
end;

CheckAndMakeDir(outJPath);
load allMis.mat

dds=8;
r=1;
nBig=768;
priorPos=[350 380 219 274 487 545
    370 215 337 545 544 344];
% priorPos=[350 380 219 274 487 545
%     370 275 337 545 544 344];
% priorPos=[350 380 219 274 487 545
%     370 220 337 545 544 344];
priorPos=[350 380 219 264 507 545
    370 170 337 555 554 344];

priorSigma=50;
% priorAmp=40;
% priorAmp=60;
nHoles=1;
locs=zeros(3000,2,5,'single');
%%
update2DFreq=5000
figure(1);
iStart=6904
for i=iStart:numel(allMis)
    mi=allMis{i};
    imgName=['Merged_sm_all/' mi.baseFilename 'ms.mrc'];
    disp([imgName '  ' num2str([nHoles i])]);
    if ~exist(imgName,'file')
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
        subplot(233);
        imags(norm);
        
    end;
    
    m8p=Crop(m8,nBig);
    % m8p=GaussHP(m8p,fHP);
    
    subplot(234);
    imags(m8);
    cc=real(ifftn(conj(fftn(ifftshift(m8p))).*fref));
    % subplot(232);
    % imags(cc);
        
    ncc=cc./(norm+k);
    
    prior=priorAmp*Gaussian(nBig,2,priorSigma,priorPos(:,iShot)');
    
    ncc=(ncc+prior).*nccMask;
    
   % subplot(235);
    % plot(sect(ncc));
    
    [mx,ix,iy]=max2d(ncc);

    subplot(232);
    imags(ncc);
    hold on;
    plot(ix,iy,'r+','markersize',10);
    hold off;
     
    model=ExtractImage(ref,[ix,iy],n8);
    subplot(235);
    
    modelx=GaussFiltDCT(model,.07)>.9; % expand by a few pixels
    
    modelImage=Downsample(modelx+60*m8,4*n8);
    jname=sprintf('%s%04d.jpg',outJPaths{iShot},nHoles);
    WriteJpeg(modelImage,jname);
    disp(jname);
    
    mi=meInsertMask(modelx,mi,2);
    WriteMiFile(mi,[outMiPath mi.baseFilename 'mi.txt']);
    
    imags(modelx+60*m8);
    
    locs(nHoles,:,iShot)=[ix iy];
    
    title(iShot);
    drawnow;
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
%     pause
end;
save(['locs' runStr '.mat'],'locs');

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
