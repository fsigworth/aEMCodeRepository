% FTEx6synth2d.m
% Example of building up a 2d image from Fourier components.
% hit space to add terms to reveal the surprise result...
% 
figure(2);
set(gcf,'menubar','none');
set(gcf,'color',[0.3 0.3 0.3]);
clf;
SetComplex;
  
n=64;
N=n*n;
p=n/2+1;
big=1e6;

sm=double(imread('smiley.tiff'));
sm=(255-sm')/255;
sz=size(sm);
ns=min(sz);  % smallest dimension of object

rm=zeros(n,n);
ctr=n/2+1;
rm(ctr-ns/2:ctr+ns/2-1,ctr-ns/2:ctr+ns/2-1)=sm(1:ns,ns:-1:1);
rm=rm-mean(mean(rm));

fm=fftshift(fftn(fftshift(rm)));
cxexp=0.5;
    pars.scl=.012;
    mulr=2000;

% subplot(2,2,1);
% imacx2(fm,cxexp,pars);

[x y]=ndgrid(-n/2:n/2-1);
r=sqrt(x.^2+y.^2);
theta=atan2(y,x)+pi;
val=-(r+theta/pi*.01);
% val(p,p)=-big;
index=0;
while index < N/8
    np=floor(index/10)+1;
    nf=ceil(100/(index+2));
    for k=1:np
        [z,i,j]=max2d(val);  % find the next point
        i2=n+2-i;
        j2=n+2-j;
        val(i,j)=-big;
        val(i2,j2)=-big;
        msk=zeros(n,'single');
        msk(i,j)=1;
        msk(i2,j2)=1;
        incr=fm.*msk;
        rIncr=real(fftshift(ifftn(fftshift(incr))));
        mysubplot(2,2,2);
        imaga(mulr*rIncr+128);
        axis off;
        pause(0.1);
        fms=fm.*(val==-big);
    mysubplot(2,2,3);
  imacx2(fms,cxexp,pars);
    title(index);
    axis off;
    pause(0.5/np);
    recon=real(fftshift(ifftn(fftshift(fms))));
    subplot(2,2,4);
    imacs(recon);
    drawnow;
    end;
    index=index+np;

end;

subplot(2,2,2);
imacs(rm);

