% DetectBadVesicles.m
fHP=.003; % Highpass, in A^-1
tRise=10;
% fHP=0;
qHP=0; % highpass leakthrough amplitude
fLP=.2;
fLP2=.4;
dispExp=.5;

cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1')
% cd('/Users/fred/EMWork/Hideki/140620')
d=dir('Info');  % entry 13 is sq10_001
id=12+9;
d(id).name
mi=ReadMiFile(['Info/' d(id).name]);
mi.basePath=AddSlash(pwd);
m=meReadMergedImage(mi);
%
n=size(m);
ds=mi.imageSize(1)/n(1);
pixA=mi.pixA*ds
figure(1);
SetGrayscale;
vm=zeros(n,'single');
% for ind=1:numel(mi.vesicle.x)
for ind=1:20
    %%
%     Pick the working (extracted) image size
    r=mi.vesicle.r(ind,1)/ds;
    nw=256;
    if r>90
        nw=384;
    end;
    pt=[mi.vesicle.x(ind) mi.vesicle.y(ind)]/ds+1;
    ipt=round(pt);
    %
    %     Subtract a single vesicle from the micrograph and display it.
    v=meMakeModelVesicles(mi,n,ind);
    figure(1);
    SetGrayscale;
    subplot(2,3,1);
     mvh=GaussHP(m-v,fHP*pixA)+qHP*(m-v);
%       mvh=GaussHP(m,fHP*pixA);
    imacs(GaussFilt(mvh,.05));
    hold on;
    plot(pt(1),pt(2),'y+');
    hold off;
    title(mi.baseFilename,'interpreter','none');
% mvh=m-v;
%
%     Extract the vesicle image
    mx=ExtractImage(mvh,ipt,nw);
    mx=GaussFilt(mx,fLP);
    subplot(2,3,2);
    imacs(mx);
    ct=nw/2+1;
    hold on;
    plot([ct ct nan 1 nw],[1 nw nan ct ct],'y-');
    hold off;
    title(['Vesicle ' num2str(ind) ';  r = ' num2str(round(r))]);
    mx=mx-median(mx(:));
    %
%     Convert to polar coordinates
    mp=ToPolar(mx);
    subplot(2,3,3);
    imacs(mp');
    title('Polar coords');
%     FT in angular space
    fmp=circshift(fft(mp'),[nw 0]);
%     Show the power spectrum
    fmp(nw+1,1:2)=0;
    smp=abs(fmp).^2;
    subplot(2,3,4);    
    imacs(-10:10,0:nw/2-1,smp(nw-9:nw+11,:).^dispExp);

    ri=round(r);
%     Determine the maximum frequency
    fw=max(2,round((r/25)));
%     Extent of radius surrounding membrane center
    mw=round((50+max(r*pixA-170,0)*.3)/pixA);
    mw=(50+max(r*pixA-170,0)*.3)/pixA;
    msk=fuzzymask(nw/2,1,mw,tRise/pixA,r);
    rmn=ri-mw+1;  % minimum radius to blank
    rmx=min(nw/2,ri+mw+1);  % max radius to blank

    hold on;
    plot(0,r,'y+');
    plot([-fw -fw fw fw -fw],[rmn rmx rmx rmn rmn],'r-');
    hold off;
    title(['Polar spectrum, dr = ' num2str(mw) ';  df = ' num2str(fw)]);
    %
    fcorr=fmp.*0;
%     fcorr(nw-fw+1:nw+fw+1,rmn:rmx)=fmp(nw-fw+1:nw+fw+1,rmn:rmx);
    fcorr(nw-fw+1:nw+fw+1,:)=fmp(nw-fw+1:nw+fw+1,:)...
        .*repmat(msk',2*fw+1,1);
%     subplot(2,3,5);
%     scorr=smp-abs(fcorr).^2;
%     imacs(-10:10,0:nw/2-1,scorr(nw-9:nw+11,:).^dispExp);
    cp=real(ifft(circshift(fcorr,[-nw 0])))';
    cr=ToRect(cp);
    subplot(2,3,5);
%     imacs(mx-cr);
    imacs(GaussFilt(mx-cr,fLP2));
    hold on
        plot([ct ct nan 1 nw],[1 nw nan ct ct],'y-');
    hold off;
    title('Corrected');
    subplot(2,3,6);
    imacs(cr);
%     plot(sect(cr),'b.-');
%     imacs(mp'-cp');
title('Residual');
    drawnow;    
    
    vm=vm+ExtractImage(cr,ipt,n,1);%%
    
end;
%%
    v=meMakeModelVesicles(mi,n);

subplot(231);
imacs(GaussFilt(m-v,.1));
subplot(232);
imacs(GaussFilt(m-v-vm,.1));
subplot(233);
imacs(vm);
