% FTEx6synth2d.m
% Example of building up a 2d image from Fourier components.
% hit space to add terms to reveal the surprise result...
% 
fs=12;


figure(2);
% set(gcf,'menubar','none');
set(gcf,'color',[0.3 0.3 0.3]);
clf;
SetComplex;
  
n=128;
N=n*n;
p=n/2+1;
big=1e6;

load TRPTop
sm=trpTop;
sm=(255-sm')/255;
sz=size(sm);
ns=min(sz);  % smallest dimension of object

rm=zeros(n,n);
ctr=n/2+1;
rm(ctr-ns/2:ctr+ns/2-1,ctr-ns/2:ctr+ns/2-1)=sm(1:ns,ns:-1:1);
rm=rm-mean(mean(rm));

fm=fftshift(fftn(fftshift(rm)));
subplot(2,2,1);
imacx(fm);

fmax=max(abs(fm(:)));

subplot(2,2,4)
imags(sm);

pause

[x y]=ndgrid(-n/2:n/2-1);
r=sqrt(x.^2+y.^2);
theta=atan2(y,x)+pi;
val=-(r+theta/pi*.01);
% val(p,p)=-big;
index=0;
while index < N/4
    np=floor(index/10)+1;
    for k=1:np
        [z,i,j]=max2d(val);  % find the next point
        i2=n+2-i;
        j2=n+2-j;
        val(i,j)=-big; 
        val(i2,j2)=-big;
        fms=fm.*(val==-big);
        fm1=zeros(n);
        fm1(i,j)=fm(i,j);
        fm1(i2,j2)=fm(i2,j2);
    end;
    index=index+np;
    subplot(2,2,3);
    imacx(fms);
    freq=round(hypot(i-ctr,j-ctr));
    title(['|f| = ' num2str(freq)],'fontsize',fs);

    term=real(fftshift(ifftn(fftshift(fm1))));
    subplot(2,2,2);
    imaga(term*1e4);
    title(num2str(fm1(i,j)/n),'fontsize',fs);
    
    recon=real(fftshift(ifftn(fftshift(fms))));
    subplot(2,2,4);
    imacs(recon);

    pause;
end;

subplot(2,2,2);
imacs(rm);

