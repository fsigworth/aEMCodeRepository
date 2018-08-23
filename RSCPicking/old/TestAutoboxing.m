miny=320;
maxy=335;
figure(4);
subplot(3,1,1);
plot(mxVar(:,miny:maxy));
axis([0 1024 0 inf]);
subplot(3,1,2);
plot(mxCC(:,miny:maxy));
axis([0 1024 0 inf]);
subplot(3,1,3);
imacs(mxCC(:,miny:maxy));
axis([0 1024 0 inf]);
figure(5);
c=mxCC(:,miny:maxy);
mxVarf=GaussFilt(mxVar,.1);
vf=mxVarf(:,miny:maxy);
v=mxVar(:,miny:maxy);
plot(c,vf,'.');
%%
figure(2);
imacs(mxVarf>200);
%%
mthresh=150;

mfc=.03;
mmin=.99;
msk0=mxVarf<mthresh;
msk=GaussFilt(msk0,mfc)>mmin;
mskCC=msk.*mxCC;

figure(2);
imacs(mskCC+.1*mxCC);
colormap jet




% Autoboxing
s0=30;  % box size
square=zeros(s0,s0);
square(1,:)=1; square(s0,:)=1;
square(:,1)=1; square(:,s0)=1;
mskr=10;
x0=100;
nx=512;
y0=0;
fc=.2;
mxc=.9;
mnc=.65;

oim=Downsample(origImg,1024);
mdois=GaussFilt(oim(x0+1:x0+nx,y0+1:y0+nx),fc);
mdo=imscale(mdois,256,.001);
mdss=GaussFilt(m(x0+1:x0+nx,y0+1:y0+nx),fc);
mds=imscale(mdss,256,.001);

cmsk=zeros(nx,nx);

cdis=mskCC(x0+1:x0+nx,y0+1:y0+nx);
cdis1=cdis;
figure(6); SetGrayscale;
di=uint8(zeros(nx,nx,3));
        cd=imscale(cdis,256,.001);
v=mxc;
nfound=0;
amps=0;
while v>mnc   && nfound<100
    [v i j]=max2d(cdis1);
    cdis1=cdis1.*(1-fuzzymask(nx,2,mskr,1,[i j]));
    point=[i j];
    if (v<mxc && v>mnc)  % accepted range of CC peaks
        nfound=nfound+1;
        amps(nfound)=v;
        cmsk=Mask(cmsk,point,1+0*square,square);
        subplot(2,2,1);
        di(:,:,1)=rot90(mdo.*(1-cmsk)+255*cmsk,1);
        di(:,:,2)=rot90(mdo.*(1-cmsk),1);
        di(:,:,3)=rot90(mdo.*(1-cmsk),1);
        image(uint8(di));
        title(nfound);
        
        subplot(2,2,2);
        di(:,:,1)=rot90(mds,1);
        di(:,:,2)=rot90(mds.*(1-cmsk),1);
        di(:,:,3)=rot90(mds.*(1-cmsk),1);
        image(uint8(di));
        title(v)

        subplot(2,2,3);
        di(:,:,1)=rot90(cd.*(1-cmsk)+255*cmsk,1);
        di(:,:,2)=rot90(cd.*(1-cmsk),1);
        di(:,:,3)=rot90(cd.*(1-cmsk),1);
        image(uint8(di));

        subplot(2,2,4);
        imacs(cdis1);
        
        drawnow;
    end;
end;
hist(amps);