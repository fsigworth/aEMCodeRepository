% rlStarV3To2.m

inStarName='particles.star';
outStarName='particles_v_V2.star';

maxNParticles=inf;

disp(['Reading ' inStarName]);
[nms,dats]=ReadStarFile(inStarName);
%%
opt=dats{1};
d=dats{2};
np=numel(d.rlnOpticsGroup);

if np>maxNParticles
    d=TrimStructFields(d,1,maxNParticles);
end;
oGroups=d.rlnOpticsGroup;
if isfield(d,'rlnOriginXAngst')
    d.rlnOriginX=d.rlnOriginXAngst/opt.rlnImagePixelSize(oGroups);
    d.rlnOriginY=d.rlnOriginYAngst/opt.rlnImagePixelSize(oGroups);
    d=rmfield(d,'rlnOriginXAngst');
    d=rmfield(d,'rlnOriginYAngst');
end;
rlnVoltage=opt.rlnVoltage(oGroups);
d.rlnSphericalAberration=opt.rlnSphericalAberration(oGroups);
d.rlnAmplitudeContrast=opt.rlnAmplitudeContrast(oGroups);
d.rlnMagnification(1:np,1)=1e4;
d.rlnDetectorPixelSize=opt.rlnImagePixelSize(oGroups);

d=rmfield(d,'rlnOpticsGroup');

%
if isfield(d,'vesicleRadius')
    d=rmfield(d,'vesicleRadius');
end;
if isfield(d,'vesiclePsi')
    d=rmfield(d,'vesiclePsi');
end;
%


disp(['Writing ' outStarName]);
WriteStarFile(nms(2),{d},outStarName);

