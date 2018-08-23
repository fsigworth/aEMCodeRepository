% reStackSimulator.m
% Note that we don't do damage/damage compensation, and don't use a
% prewhitening filter.

sim=struct;
sim.n=128;

sim.nmi=64; % no. of micrographs
sim.npm=64;  % particles per micrograph
% sim.nmi=64; % no. of micrographs
% sim.npm=16;  % particles per micrograph
sim.mapName='KvMap.mat';
def=[2.5 4]; % defocus values: min, max, offset if 2nd exp
weights=1;
doses=50;
sim.mergeMode=3;

sim.pIso=.3;
sim.angleLimits=[0 20 90];
sim.mbnOffsetA=48;
sim.ds=2;
sim.sigmaC=0;
sim.sineWeighting=1;

sim.sigmaSVesicle=.1;
sim.rVes=[150 300];  % min, max (uniform)
sim.sVes=[.003 .01]; % min,max (uniform)
sim.aParticle=[.2 .3];
sim.usePWFilter=1;

mi0=reSimMiStructure(weights*0,doses);
mi0.weights=weights;
mis=cell(sim.nmi,1);
for i=1:sim.nmi
    d1=def(1)+(def(2)-def(1))*rand;  % random defocus values
    mi0.ctf(1).defocus=d1;
    mi0.ctf(1).B=40;
%     mi0.ctf(2).defocus=d1+def(3);
%     mi0.ctf(2).B=80;
    mi0.weights=weights;
    mis{i}=mi0;
end;

[si,imgs0]=reMakeFakeImages(sim,mis,sim.mapName);
%%
nim=sim.nmi*sim.npm;
sim.sigmaN=.14;

% Store the simulation parameters
si.sim=sim;
n=sim.n;
noise=sim.sigmaN*randn(n,n,nim);
imgs=imgs0+noise;

figure(10);
ImagicDisplay3(BinImage(Crop(imgs,96,1),4));

name0=sprintf('Sim%dn%03dk%dm%dtsi.mat',n,round(100*sim.sigmaN),round(nim/1024),sim.mergeMode);
% name0=['Sim' num2str(sim.n) 'n007k8tsi.mat'];
disp('Getting a stack file name');
[name,pa]=uiputfile(name0);
if isnumeric(name)
    return
end;
cd(pa);
name=[pa name];
disp(['Writing ' name]);
save(name,'si');
stackName=[name(1:end-6) 'stack.mrc'];
disp(['Writing ' stackName]);
WriteMRC(imgs,si.pixA,stackName);
disp('done.');
