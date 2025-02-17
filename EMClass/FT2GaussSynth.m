% FT2GaussSynth.m
% Example of building up a 2d image from Fourier components.
% 
vlSetDarkGraphics(18,[],1);

figure(1);

clf;
nx=100;
nx2=2*nx+1;
ndis=25;
ctr=nx+1;
d=.5;
a1=2;
a2=2;
xLim=1;
[x,y]=ndgrid(-xLim:1/nx:xLim);  % 201 points
m=a1*a2*(exp(-pi*((a1*x).^2+(a2*(y+d)).^2))+exp(-pi*((a1*x).^2+(a2*(y-d)).^2)));
xdis=x(ctr-ndis:ctr+ndis);
subplot(231);
imags(x(:,1),y(1,:),m);


fm=fftshift(fftn(ifftshift(m)))/nx^2;

fm=real(fm);

subplot(232);
pars.scl=10;
imags(xdis,xdis,Crop(fm,ndis));

fm1=fm;
nPoints=20;
termSum=zeros(nx2,nx2,'single');
fTermSum=zeros(nx2,nx2,'single');
for i=1:nPoints
    [mx,i1,j1]=max2d(abs(fm1));
    i2=2*ctr-i1;
    j2=2*ctr-j1;
    fTerm=zeros(nx2,nx2,'single');
    val=fm(i1,j1);
    fm1(i1,j1)=0;
    fm1(i2,j2)=0;
    fTerm(i1,j1)=val/mx;
    fTerm(i2,j2)=conj(val)/mx;
    term1=real(fftshift(ifftn(ifftshift(fTerm))))*1e4;
    subplot(234);
    imaga(mx*term1*128/nx2^2+128);
    
    termSum=termSum+mx*term1;
    fTermSum=fTermSum+fTerm*mx;
    subplot(235);
    imaga(Crop(fTermSum,ndis)*128+128);
    
    subplot(236);
    imags(termSum);
    pause;
end;
% subplot(223);
% plot(x(:,1),sect(fm'));
% 
% subplot(224);
return

mx=Crop(1-m.smiley1k,2048);
m=DownsampleGeneral(mx,n,0.2);
m(m>1.13)=1.13;
mf=GaussFilt(m,.3);

imags(mf);
% plot(sect(mf))


fm=fftshift(fftn(fftshift(mf)));
fms=fm.*0;
subplot(2,2,1);
pars.scl=.001;
pars.sat=.7;
scl=imacx2(Crop(fm,ndis),1,pars);

% [x y]=ndgrid(-n/2:n/2-1);
% r=sqrt(x.^2+y.^2);
% theta=atan2(y,x)+pi;
[r,theta]=Radius(n);
val=-(r+theta/pi*.01);
active=false(n);
% val(p,p)=-big;
index=0;
while index < n/4
    np=floor(index/10)+1;
    for k=1:np
        [z,i,j]=max2d(val);  % find the next point
        i2=n+2-i;
        j2=n+2-j;
        mft=zeros(n);
        mft(i,j)=mf(i,j);
        mft(i2,j2)=mf(i2,j2);
        val(i,j)=-big;
        val(i2,j2)=-big;
        mfs=mfs+mft;
    end;
    subplot(222);
    term=real(fftshift(ifftn(ifftshift(mft))));
    imaga(256*term);
    title(mf(i2,j2));
    index=index+np;
    subplot(2,2,3);
  imacx2(Crop(mfs,ndis),1,pars);
    title(index);

    recon=real(fftshift(ifftn(fftshift(mfs))));
    subplot(2,2,4);
    imaga(recon*255);

    pause;
end;

subplot(2,2,2);
imags(mf);

