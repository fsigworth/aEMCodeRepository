% CTFTestMicrograph.m
% Make a fake micrograph for checking the ctffind3 parameters




defs=[2 2 2];
deltadef=.5;
thetas=[0 pi/10 pi/4];
n=4096;
Cs=2;
B=100;
alpha=.02;
sigmaN=1;
pixA=1.5;

m0=randn(4096);
m1=randn(4096);
fm=fftn(m0);

for i=1:numel(defs)
%%
%     i=2;
%     thetas(2)=30*pi/180;
%     deltadef=.5;
ct=CTF(n,pixA,EWavelength(200),defs(i),Cs,B,alpha,deltadef,thetas(i));
m=real(ifftn(fm.*ifftshift(ct)))+sigmaN*m1;

sp=RadialPowerSpectrum(m,0,4);
plot(sp);

name=['TestMicrograph' num2str(i) '.mrc']
WriteMRC(m,pixA,name);
%
P.lambda=EWavelength(200);
P.defocus=1:.2:3;
P.deltadef=-1:.2:1;
P.theta=0:.2:pi/4;
P.alpha=alpha;
P.B=B;
P.Cs=Cs;

[Ps, c, sp]=FitCTF(m,P,pixA,8,100,1);
opts.blockSize=512;
Ps2=FitCTFk2(m,P,pixA,1,[.02 .12],opts);
Ps.theta=Ps.theta*180/pi;
Ps2.theta=Ps2.theta*180/pi;

Ps
Ps2
end;

% CTFFIND3 results
% 25000 15000 0
% 24000 15000 18.3
% 24000 15000 45