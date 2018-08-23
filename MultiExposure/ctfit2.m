function [P, c]=ctfit2(mc,Pa,pixA,maxres,minres,test)
% function [P c]=ctfit2(mc,Pa,pixA,maxres,minres)
%  Derived from ppcfit2, uses the algorithm of N. Grigorieff's CTFFIND3.
% Determine the CTF parameters P from the image mc.  The cross-correlation
% coefficient c is returned too.
% Positive defocus values correspond to *underfocus*.
% Parameter ranges to be searched are defined by haveing ranges of values in
% the input structure Pa.  For example,
% Pa.lambda=EWavelength(200);
% Pa.defocus=0:0.5:20;
% Pa.deltadef=-2:.2:2;
% Pa.theta=0;
% Pa.alpha=0:.03:.12;
% Pa.Cs=2;
% Pa.B=1000;
% maxres and minres are in A.
% FinalPars=ctfit2(mc,Pa,1.4,8,50);
% The test argument controls display of results.  test=1 means graphical
% output (default).  test=2 gives text display of execution times.

disp('ctfit2');

if nargin<6
    test=0;
end;
NumRefs=3;
persistent R;
overfocus=0;

mc=mc-mean(mc(:));
n=size(mc,1);
%
% maxres=6;
% % maxres=10;
% minres=40;
% minres
res=pixA; % pixel size in A.
Pa.res=res;

% nu=128;  %  FFT block size, must be a multiple of 4.
nu=n/4;
nu=n/2;

w0=100;  % effective res of filter, in A.
disexp=0.2; % display exponent
pwExp1=4;  % pre-whitening positive exponent
pwExp4=2; % pre-whitening negative exponent (r^4)

w=nu*res/w0;  % filter radius.  Niko has it = (1/40A) in freq. space.


df=1/(res*nu);  % spatial frequency unit
freqs=(0:df:df*(nu/2-1))';

% Compute local power spectra
[nx ny]=size(mc);
nv=round(nu/4);   % overlap
tx=TileCoords(nx,nu,nv);  % no. of tiles in X direction
ty=TileCoords(ny,nu,nv);  % no. of tiles in Y direction

if test
    figure(1);
    SetGrayscale;
    subplot(2,3,1);
    imacs(mc);
    LabelImage(mc,res/10,'nm');
    drawnow;
end;

sd=zeros(tx,ty);
% Make a window for spectral estiation
window=fuzzymask(nu,2,0.45*nu,0.1*nu);
winsum=sum(window(:));

sps=zeros(nu,nu,tx,ty);

if test>1
    disp('Computing spectra...');
    tic
end;
for ix=1:tx
    [x0 x1 u1]=TileCoords(nx,nu,nv,ix);
    
    for iy=1:ty
        [y0 y1 v1]=TileCoords(ny,nu,nv,iy);
        % I don't think it's necessary to remove gradients, but Niko does
        % this:
        tm=RemoveGradients(double(mc(x0+1:x0+nu,y0+1:y0+nu)));
        tm=tm.*window;
        tm=tm-sum(tm(:))*window/winsum;
        
        % subplot(2,1,1); imacs(tm); title([ix iy]);
        
        sp2=abs(fftn(tm)).^2;
        % subplot(2,1,2); imacs(fftshift(sqrt(sp2)));
        % drawnow;
        sps(:,:,ix,iy)=sp2;
        sd(ix,iy)=sqrt(sum(sp2(:)))/nu;
    end;
end;
if test>1
    toc
end;
% subplot(2,3,1);
% sumsp=sum(sum(sps,4),3);
% imacs(fftshift(sumsp.^disexp));
% LabelImage(nu,df,'A^{-1}');


% Make a histogram of the s.d. values and find tiles having s.d. near the
% mode.
nbins=tx+2;
hthresh=0.2;
% hthresh=-1;  % no selection

sdmn=min(sd(:));
sdmx=max(sd(:));
dbin=(sdmx-sdmn)/(nbins-1);
bins=sdmn:dbin:sdmx;

[h x]=histc(sd(:),bins);

[mxh j]=max(h);
% j is the mode of the distribution
% Find the first bin below hthresh
jmax=j-1+find(h(j:nbins)<hthresh*mxh,1,'first');
if numel(jmax)<1
    jmax=nbins;
end;
sdmax=x(jmax);
% find the first bin below hthresh
jmin=find(h(1:j)<hthresh*mxh,1,'last');
if numel(jmin)<1
    jmin=1;
end;
sdmin=bins(jmin);
sdmax=bins(jmax)+dbin;

% Show histogram
% subplot(2,3,2);
% bar(x,h);
% hold on
% plot(sdmin,0,'w.');
% plot(sdmax,0,'w.');
% hold off

% Show regions that we used.
subplot(2,3,2);
imacs((sd<=sdmax).*(sd>=sdmin));
drawnow;

cumsp=zeros(nu,nu);
count=0;
for ix=1:tx
    for iy=1:ty
        if (sd(ix,iy)<=sdmax) && (sd(ix,iy)>=sdmin)
            cumsp=cumsp+sps(:,:,ix,iy);
            count=count+1;
        end;
    end;
end;
if test>1
    disp([num2str(count) ' tiles of ' num2str(tx*ty) ' used']);
end;
% pre-whitening correction
r=fftshift(Radius(nu))/nu;
cumsp0=cumsp; %%%
cumsp=cumsp.*exp(r.*(pwExp1-pwExp4*r.^3));

% remove the crystal spots
% cumsp=fftshift(RemoveSpots(fftshift(cumsp),[3.5 3]));

% Filter the square root of the spectrum.
sqrsp=sqrt(cumsp);
subplot(2,3,4);
imacs(fftshift(sqrsp.^(disexp*2)));
% % LabelImage(nu,df,'A^{-1}');

kernel=fuzzymask(nu,2,w,w/10);
kernel=kernel/sum(kernel(:));
% convolve with kernel
filtsp=real(ifftn(fftn(sqrsp).*fftn(fftshift(kernel))));

diffsp=fftshift(cumsp-filtsp.^2);
% diffsp=fftshift(cumsp./filtsp-filtsp);
% diffsp=fftshift(cumsp);

% look at the subtracted spectrum

% % subplot(2,3,3);  % Compare the smoothed and original spectra
% % radSpecs=[Radial(fftshift(filtsp))' Radial(fftshift(sqrt(cumsp)))'];
% % semilogy(freqs,radSpecs);
% % xlabel('A^{-1}');
% % drawnow;
disp('Radial');

radsp=Radial(fftshift(cumsp));


% determine resolution limits

rmin=1/(df*minres);
rmax=1/(df*maxres);
outerlimits=fuzzymask(nu,2,rmax,rmax/10);
limits=outerlimits-fuzzymask(nu,2,rmin,1);

subplot(2,3,5);
imacs(limits.*diffsp);
% % LabelImage(nu,df,'A^{-1}');

subplot(2,3,6);
rdiff=Radial(diffsp.*limits);
% % plot(freqs,rdiff);
drawnow;


% Set up for computing references
% --search limits defined here.
% Pa.lambda=EWavelength(300);
% Pa.defocus=0:0.2:6;
% Pa.deltadef=-.6:.2:.6;
% Pa.theta=0:pi/8:3*pi/8;
% Pa.alpha=0:pi/8:pi;
% Pa.Cs=4;
% Pa.B=100;
% Pa.res=res;

Pao=Pa;
Pao.defocus=-Pa.defocus;
% Pa.B=200;

halflimits=limits(nu/2+1:nu,:);
halflimits(1,nu/2:nu)=0;  % zero out the upper meridian

halfsp=diffsp(nu/2+1:nu,:).*halflimits;
halfspv=halfsp(:);
halfspv=halfspv/sqrt(halfspv'*halfspv);

% % Here we test for an existing set of references, and compute new ones
% % only if necessary.
% First, check whether Pa0 exists at all.
b=1;
try
    b=numel(R);
catch
    b=0;
end;
% Now see if the problem has changed.  When we change this script to a
% function, we can use persistent variables for refs, refsz, Pa0.
ind=0;
for i=1:b
    if StructsEqual(R(i).Pa0,Pa)
        ind=i;
    end;
end;

if ind==0
    ind=min(b+1,NumRefs);
    disp('Making new CTF references');
    [R(ind).refs R(ind).refsz]=MakeCTFRefs(nu,res,Pa,halflimits);
    [R(ind).refso R(ind).refsz]=MakeCTFRefs(nu,res,Pao,halflimits);
    
    R(ind).refs=reshape(R(ind).refs,nu^2/2,prod(R(ind).refsz));
    R(ind).refso=reshape(R(ind).refso,nu^2/2,prod(R(ind).refsz));
    R(ind).Pa0=Pa;
%     whos R
end;

if test>1
    disp('Cross-correlation')
end;

if overfocus
    cc=halfspv'*R(ind).refso;
else
    cc=halfspv'*R(ind).refs;  % Cross correlation done here!
end;
sz=size(cc);
[mxc mxi]=max(cc);
% mxi
[m l i j k]=ind2sub(R(ind).refsz,mxi);

if test>1
    Initial_c=mxc
end;

if overfocus
    Ps=Pao;
    P=Ps;
else
    Ps=Pa;
    P=Ps;
end;

P.defocus=Ps.defocus(i);
P.deltadef=Ps.deltadef(j);
P.theta=Ps.theta(k);
P.alpha=Ps.alpha(l);
P.B=P.B(m);
% P

% Make the split-power spectrum display
subplot(2,3,5);
dspect=imscale(limits.*diffsp);
dref=imscale(limits.*CTF(nu,res,P).^2);
dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
imac(dspect);
LabelImage(nu,df,'A^{-1}');

% % Make the radial display
% subplot(2,3,6);
% plot(freqs,[Radial(imscale(limits.*diffsp))' Radial(dref)']);
% legend('data','fit');
% xlabel('A^{-1}');

drawnow;

ccx=reshape(cc,R(ind).refsz);


% figure(3); clf;
% colormap jet
% make a contour plot of alpha vs defocus
subplot(2,3,3);
if numel(Ps.alpha)>1
    contourf(Ps.alpha,-Ps.defocus,squeeze(ccx(m,:,:,j,k))',10);
    xlabel('Alpha');
    ylabel('Defocus');
    colorbar;
    title('Underfocus is negative');
end;
subplot(2,3,3);
contourf(Ps.deltadef,-Ps.defocus,squeeze(ccx(m,l,:,:,k)),10);
xlabel('Delta-defocus');
ylabel('Defocus');
colorbar;
drawnow;

% Do the optimization
% disp('Final optimization');
p=[P.defocus; P.deltadef; P.theta];

if test>1
    disp('Final optimization.');
    disp('Initial values:');
    disp([-p(1) p(2)]);
end;

P.alpha=.07;
P.theta=0;
for is=1:2  % do extra restarts of the optimization, to get out local optima.
    %     disp(p');
    p=Simplex('init',p,[0.1 0.1]);
    for i=1:80
        P.defocus=p(1);
        P.deltadef=p(2);
        c=CTFCrossCorr(diffsp,res,limits,P,0);
        p=Simplex(-c);
        %               if mod(i,10)==0
        %                 p'
        %            end;
    end;
    % disp([p(1:2)' 180/pi*p(3:4)']);
end;

subplot(2,3,4);
% Old circularly-averaged display
% dspect=imscale(fftshift(sqrsp.^(disexp*2)));
% dspect=(imscale(ImEqualize(GaussFilt(fftshift(sqrsp),0.2)).^2)*1+0);
% [dref mr ma]=imscale(limits.*CTF(nu,res,P).^2);
% dref=ma+mr*CTF(nu,res,P).^2;
% % dref=ImEqualize(CTF(nu,res,P).^2);
% dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
% imac(dspect);
% LabelImage(nu,df,'A^{-1}');
dspect=imscale(limits.*diffsp);
dref=imscale(limits.*CTF(nu,res,P).^2);
dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
xs=-nu/2:nu/2+1;
imac(xs*df,xs*df,dspect);
xlabel('A^{-1}');
title(['\Deltaz=' num2str(-p(1)) '   Astig=' num2str(p(2))...
    '   \theta=' num2str(360/pi*p(3))]);
title(['CC =' num2str(c)]);

subplot(2,3,5);
dspect=imscale(limits.*diffsp);
dref=imscale(limits.*CTF(nu,res,P).^2);
dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
imac(dspect);
LabelImage(nu,df,'A^{-1}');
title(['\Deltaz=' num2str(-p(1)) '   Astig=' num2str(p(2))...
    '   \theta=' num2str(360/pi*p(3))]);

% Make the radial display
subplot(2,3,6);
rs=Radial(limits.*diffsp);
rs=150*(rs-min(rs))/(max(rs)-min(rs));
radCorrSpecs=[rs Radial(dref)];
plot(freqs,radCorrSpecs);
xlabel('A^{-1}');
drawnow;

% P.defocus=-P.defocus;  % change the polarity to match scope.
if test>1
    c
    Final_values=[P.defocus P.deltadef 180/pi*P.theta P.alpha]
end;
% CTF_Parameters=P


