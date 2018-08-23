% rlTestUnblurFilter

cd('/Users/fred/EMWork/Hideki/170810-vpp');
runDir='MC2Run2/';
d=dir(runDir);
n=[3840 3840];
% annulus=1-(fuzzymask(n,2,1870,0)-fuzzymask(n,2,1910,0));


mTemp=ReadMRC([runDir d(3).name]);  % summed stack
% sz=size(mTemp);
% annulus=1-(fuzzymask(sz,2,0.9*sz,0)-fuzzymask(sz,2,0.98*sz,0));
% q=annulus.*mTemp;
% me=sum(q(:)/sum(annulus(:));
m0=Crop(mTemp,n,0,mean(mTemp(:)));  % summed stack
mTemp=ReadMRC([runDir d(4).name]);  % dose-weighted stack
md=Crop(mTemp,n,0,mean(mTemp(:)));
disp('Reading movie');
%%
d1=dir;
mTemp=single(ReadMovie(d1(4).name));

disp('Downsampling');
mv=Downsample(Crop(mTemp,2*n,1,mean(mTemp(:))),n,1);

%
disp('Spectrum');
sp0=RadialPowerSpectrum(m0);
spd=RadialPowerSpectrum(md);
mvs=RemoveOutliers(sum(mv,3));
sps=RadialPowerSpectrum(mvs); % raw stack
%
% imags(annulus);
%
disp('Filter');
H=rlUnblurFilter([n 24],2,2.09,300);
%%
sp0c=sp0;
sp0c(sp0>5)=5;
spdc=spd;
spdc(spd*25>5)=5/25;
spsc=sps;
spsc(sps>5)=5;
figure(1);
plot([sp0c spdc*25 spsc sectr(H).^2/10]);
legend({'unweighted' 'dose-weighted' 'stack sum' 'H^2'});
