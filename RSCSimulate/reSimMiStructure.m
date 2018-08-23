
function mi=reSimMiStructure(defoci,doses)
% Make a simple mi structure for simulations
if nargin<1
    defoci=[0 0];
end;
if nargin<2
    doses=[54 30];
end;
nCTFs=numel(defoci);
mi=meCreateMicrographInfoStruct14;
mi.imageSize=[3840 3840];
mi.pixA=1.247;
mi.kV=200;
mi.cpe=0.8;
mi.camera=5;
mi.frameDose=3*ones(1,30);
mi.frameSets=[2 19; 21 30];
mi.doses=doses;
mi.weights=ones(1,nCTFs);
mi.ctf=struct;

mi.ctf.lambda=.0251;
mi.ctf.defocus=defoci(1);
mi.ctf.deltadef=0;
mi.ctf.theta=0;
mi.ctf.alpha=.02;
mi.ctf.Cs=2;
mi.ctf.B=40;
mi.ctf.ampFactor=1;
for i=2:numel(defoci)
    mi.ctf(i)=mi.ctf;
    mi.ctf(i).defocus=defoci(i);
end;
mi.damageModelCode='.184*f.^-1.665+2.1;';
mi.noiseModelPars=[1.9612 0 0 1 100 0.0535 2 1];
mi.noiseModelCode={'ag=p(1);'
    'af1=p(1);'
    'af2=p(2);'
    'ag=p(3);'
    'sigma=p(4);'
    'bf=p(5);'
    's0=p(6);'
    'f1exp=p(7);'
    'f2exp=p(8);'
    'gauss=exp(-f.^2/(2*sigma^2));'
    'f1=(.01./f).^f1exp;'
    'f1(1)=f1(2);'
    'f2=(.01./f).^f2exp;'
    'f2(1)=f2(2);'
    'b2=1./((f/bf).^2+1);'
    'spec=max(0,(ag*gauss+af1).*f1+af2*f2);'
    'shot=max(0,s0*b2);'};



