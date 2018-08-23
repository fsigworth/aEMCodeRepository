% TestAnalyzeK2Movies

cd('/Volumes/TetraSSD/130822/GluR2GFP_glutamate_slot4/movie_frames/');
movName='Aug23_00.09.04.mrc';
movName='Aug22_22.31.49.mrc';
movName='Aug22_21.03.26.mrc';
movName='Aug22_19.51.13.mrc';
movName='Aug23_02.32.07.mrc';
r1Name='CountRef_Aug23_00.09.04.dm4';
r2Name='CountRef_Aug23_02.02.35.dm4';

disp(['Reading ' movName]);
mv=ReadMRC(movName);
disp(['Reading ' r1Name]);
r1=ReadEMFile(r1Name);
disp(['Reading ' r2Name]);
r2=ReadEMFile(r2Name);

%%
% See which way the reference works
nr=1;  % this is best
m0=sum(single(mv),3);
% m0=m0-mean(m0(:));
rr=rot90(r1,nr);
p=m0(:)'*rr(:)
%%


%%
figure(1);
SetGrayscale;
mf=single(mv);
xs=0:100;
h=hist(mf(:),xs);
bar(xs,h);
drawnow;
hs=sum(h);
%%
cutoff=find(h<hs*1e-7,1,'first')
mf(mf>=cutoff)=0;
mfs=sum(mf,3);

%%
nim=size(mf,3)
mf0=mf-mean(mf(:));
mfs0=sum(mf0,3);
% acf=fftshift(abs(fftn(double(mfs0))).^2);
% aci=fftshift(sum(abs(fft2(double(mf0))).^2,3));
% acd=acf-nim*aci;
return

%%
nim=size(mf,3);
n=size(mfs0);
csp0=zeros(n);
 for i=1:nim
% i=12;
csp0=csp0+fftn(mf0(:,:,i)).*conj(fftn(mfs0-mf0(:,:,i)));
% end;

 msk=ones(n);
msk(1,120:120:n(2))=0;
msk(1,119:120:n(2))=0;
msk(1,1:120:n(2))=0;
csp=csp0.*msk;


acs=fftshift((ifftn(csp)));
cacs=Crop(acs,256);
% cacs(129,129)=cacs(1);
% imacs(cacs);
%
cacc=cacs;

%  Hardwired peak correction
us=-.21;
us2=-.08;
mul=[1 1 1 1 1; 1 0 0 0 1; 1 1 0 1 1; 1 0 0 0 1; 1 1 1 1 1]';
pmul=Crop(mul,256);
frac=[us us2 us; 0 1 0; us us2 us]';
pfrac=Crop(frac,256);
pkv=cacc(129,129)-(pmul(:)'*cacc(:))/sum(mul(:));
cacc=cacs-pkv*pfrac;
imacs(Crop(cacc,64));
imacs((cacc));
title(i);
drawnow;
 end;
 return
 
 % q=Crop(GaussFilt(cacc,.1),64);
% contourf(q',40);
% plot(sect(q))
% imacs(q);
%  imacs(Crop(GaussFilt(cacc,.1),64));
%%
pk=cacs;
ct=129;
t=2;
pk(ct-t:ct+t,ct-t:ct+t)=mean(cacs(:));
for i=1:40
    fpk=GaussFilt(pk,.1);
    pk(ct-t:ct+t,ct-t:ct+t)=fpk(ct-t:ct+t,ct-t:ct+t);   
    imacs(pk);
    title(i);
    pause(.1);
end;



%%
rcsp=real(csp);
q=fftshift(rcsp);
qmax=1e11;
q(q>qmax)=qmax;
imacs(q)

plot(sect(q'));
r=ifftshift(q');
rf=GaussFilt(r,.01);
plot(rf(:,1));
% plot(sect(q'));


return

%%
figure(2)
plot(sect(cacs));

figure(1);
csps=fftshift(csp);
bcsps=BinImage(real(csps),4);
rcsps=Radial(bcsps);
semilogy(rcsps);
% 
% plot([sect(GaussFilt(real(csps),.01)) sect(GaussFilt(imag(csps),.01))]);