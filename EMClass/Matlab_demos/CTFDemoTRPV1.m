% CTFDemoTRPV1.m
rotAngle=45;

[m,s]=ReadMRC('emd_5778.map');
if mod(rotAngle,90)~=0
    disp('Rotating...');
    mr=rsRotateImage(m,rotAngle);
    disp('done.');
else
    mr=m;
end;
%%
m1=squeeze(sum(mr,1)); % projection
n=size(m1,1);
kV=300;
def=2;
Cs=2;
B=30/4;
% B=0;
alpha=.02;
deltadef=0;
theta=pi/4;
c=CTF(n,s.pixA,EWavelength(kV),def,Cs,B,alpha,deltadef,theta);
ndis=n/2;
mul=1/1.77;
% mul=.2;
add=130;
% sigma=0;

% set(gcf,'color',.3*[1 1 1]);
% vlSet1080Figure
vlSetDarkGraphics
figure(1);

% SetComplex;
mysubplot(2,3,1);
imags(m1);
text(0,0,'Projection',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;
mysubplot(2,3,4);
fm1=fftshift(fftn(m1));
imacx2(Crop(fm1,ndis),.4);
text(0,0,'FT of Projection',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;
mysubplot(2,3,5);
imacx2(Crop(c,ndis));
text(0,0,'CTF',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;
mysubplot(2,3,2);
imags(fftshift(real(ifftn(ifftshift(c)))));
text(0,0,'Point-spread function',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;
mysubplot(2,3,6);
imacx2(Crop(fm1.*c,ndis),.4);
text(0,0,'FT of image',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;
mysubplot(2,3,3);
mc=real(ifftn(ifftshift(fm1.*c)));
imaga(mc/1.77+128);
text(0,0,'Image',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;
% imags(mc);
return
%%
sigma=30;

defs=[1 2 3];
nc=numel(defs)+1;

figure(2);
% set(gcf,'color',.3*[1 1 1]);
subplot(3,nc,1);
szs=(1:n)*s.pixA;
imags(szs,szs,m1);
xlabel('angstroms');
subplot(3,nc,nc+1);
imags(szs,szs,-m1);
xlabel('angstroms');
for j=2:nc
    ind=j-1;
    N=randn(n)*sigma;
    N1=fftn(N);
    def=defs(j-1);
    c=CTF(n,s.pixA,EWavelength(kV),def,Cs,B,alpha,deltadef,theta);
    mc=real(ifftn(ifftshift(fm1.*c)));
    subplot(3,nc,j);
    imaga((mc+N)*mul*.7+add);
    axis off;
    title([num2str(def) '\mum']);
    mcf=real(ifftn(ifftshift(-fm1.*abs(c))));
    subplot(3,nc,j+nc);
    imaga((mcf+N)*mul*.7+add);
    axis off;
    if j==2
        title('Phase-flipped');
    end;
    k=(sigma^2*2e-5)+1e-5;
    wf=(fm1.*abs(c)+N1).*abs(c)./(k+abs(c).^2);
    mwf=-real(ifftn(ifftshift(wf)));
    subplot(3,nc,j+2*nc);
    imaga(mwf*mul/1.5+add+100/(sigma+2));
    axis off;
    if j==2
        title('Wiener filtered');
        cs=c;
        mcs=mc;
    else
        cs(:,:,ind)=c; % accumulate the ctfs
        mcs(:,:,ind)=mc;
    end;
    
end;

% Compute a combined Wiener-filtered image


