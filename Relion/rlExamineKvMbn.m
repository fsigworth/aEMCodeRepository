% ExamineKvMbn.m

ds=1;

cd('/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/')
msk=ReadMRC('Stack2/KvMask192.mrc');
n=size(msk,1);
mskd=Downsample(msk,n/ds);
mskr=rsRotateImage(mskd,-18);

% cd('PostProcess/job013')
cd('Stack3/PostProcess/job023')

[m,s]=ReadMRC('postprocess.mrc');
n=size(m,1);
pixA=s.pixA*ds;
ns=n/ds;
md=Downsample(m,ns);
mdr=rsRotateImage(md,-18);
figure(1);
% ShowSections(mdr,[48 68 65],45)
ShowSections(mdr,[],45)

drawnow;
%% 
figure(1);
offset=2/ds*[48 68 65];
ShowSections(mdr,offset,45);
figure(2);
ShowSections(mdr.*(1-mskr),2/ds*[48 68 65],45);
%%
disc=fuzzymask(ns,2,55/ds,2);
cyl=repmat(disc,1,1,ns);
sph=fuzzymask(ns,3,110/pixA,2);
cmsk=(cyl-mskr).*sph;

ShowSections(cmsk.*mdr,offset,45)
% ShowSections(cmsk,offset,45);
%% Get the expectation along z

mzSum=squeeze(sum(sum(cmsk)));
mbSum=squeeze(sum(sum(cmsk.*mdr)));

figure(3);
clf;
plot([mzSum 1000*mbSum]);

mzSum=max(mzSum,10);

plot(mbSum./mzSum);


%%
% Create a filter
a0=.02;
a1=.1;
d=1;
lambda=.025;
f=Radius3(ns)/(pixA*n);
r=pi * d * 1e4 * lambda * f.^2;  % r=pi when df=.06

[X,Y,Z]=ndgrid(-ns/2:ns/2-1);
kZ=0;
kZ=.001;
z0=40/ds;
z1=75/ds;
zRamp=max(0,Z-z0)*kZ;
zRamp(Z>z1)=0;

% function is the ratio  (.1 + \pi d f^2) / (.02 + \pi d f^2)

cT=(a1+r)./(a0+r);  % CTF inverse filter

hT=1;
% f1=f;
% f1(f1==0)=1;
% fh=.0002*ds;
% k=1;
% hT=exp(k*(fh./f1).^2);

mdf=real(ifftn(fftn(mdr).*ifftshift(hT.*cT)))+zRamp;
figure(2);
ShowSections(mdf,[],45);

WriteMRC(mdf,pixA,'postprocessIFilt1.mrc');

mzSum=squeeze(sum(sum(cmsk)));
mbSum=squeeze(sum(sum(cmsk.*mdf)));

figure(3);
clf;
plot([mzSum 1000*mbSum]);

mzSum=max(mzSum,10);

plot(mbSum./mzSum);

%%
figure(4);
imagsar(mdf);
figure(5);
imagsar(mdf.*mskr);