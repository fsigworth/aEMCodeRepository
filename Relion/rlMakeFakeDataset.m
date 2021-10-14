% rlMakeFakeDataset
% Comparing our angle assignments with Relion's
% Here is the result:
% rl_rot = -phi = -angs(i,1)
% rl_tilt=theta= angs(i,2)
% rl_psi = -psi-90 = -angs(i,3)-90;
% we apply shifts after rotating and projecting.
% The star file written out here corresponds to the data.star file produced in
% 3DRefine.

% We make nAmps cycles through all the angles, but with different amplitude
% values given by the amps vector.


% mapName='/Users/fred/aEMCodeRepository/AMPAR/KvMap.mat';
disp(['Working directory: ' pwd]);
p=struct;
p.baseName='SimTMStk1uDef04Ampa';
% p.baseName='SimStack1uDef01Ampd';
p.amps=[.04];

pa=fileparts(which('arGetRefVolumes'));
p.mapName=[pa '/KvMap.mat'];
p.dotCount=100; % number of computed templates per dot printed out.

% p.outDir='/Users/fred/EMWork/Simulations/Relion/';
% p.outDir='/Users/fred/EMWork/Simulations/Relion2/Fred2_HRPicking/SimStacks';

% p.outDir='/Users/fred/EMWork/Simulations/Relion/';
p.outDir='SimTMStacks/';

stackName=[p.baseName '.mrcs'];
starName=[p.baseName '.star'];
logName=[p.baseName 'Log.txt'];
refName=[p.baseName 'Ref.mrc'];
nAmps=numel(p.amps);
p.defMin=1;  % Assign min, max defocus
p.defMax=1.5;
p.imgsPerMicrograph=200;
p.kV=300;

% p.imgSize=144; % for full Kv complex
% p.mapClip=0; % blank the map from z=1:mapClip
% p.mapZShift=0;

% To make isolated TM region
p.imgSize=128;
p.mapClip=74; % blank the map from z=1:mapClip
p.mapZShift=-24;

ds=1;
p.useWhiteNoise=1; % Lorentzian noise
p.sigma=11.925; % noise makes unity variance (empirical)
p.sigma=1;
p.useUniqueNoise=1;

p.B=60;       % ctf B factor

disp('MakeFakeDataset:');


s=load(p.mapName);  % gets s.map, s.pixA; map is 108^3 in size.
% ShowSections(s.map);

m=DownsampleGeneral(s.map,p.imgSize,1/ds);
p.pixA=s.pixA*ds;
if p.mapClip>0
    msk=ones(p.imgSize,p.imgSize,p.imgSize,'single');
    msk(:,:,1:p.mapClip)=0;
    fmsk=GaussFilt(msk,.1);
    m=circshift(m.*fmsk,[0 0 p.mapZShift]);
end;


p.symmetry=4;
p.angStep=2;
%     angStep=6; %%%
p.psiStep=90; % all angles are in degrees
p.shiftMag=1;

p.nTheta=round(180/p.angStep)+1;
dTheta=180/(p.nTheta-1);
p.nPsi=ceil(360/p.psiStep);
dPsi=360/p.nPsi;
psis=(0:dPsi:360-dPsi)';
p.maxNPhi=ceil(360/(p.angStep*p.symmetry))';

nAngs=0;
angs=zeros(0,3);
for i=1:p.nTheta
    theta=mod(90+(i-1)*dTheta,180); % start at 90 and go up
    nPhi=ceil(sind(theta)*p.maxNPhi); % Sample phi sparsely when sin(theta) is small.
    dPhi=360/(nPhi*p.symmetry);
    for j=1:nPhi
        phi=(j-1)*dPhi;
        % enter all the psi values
        angs(nAngs+1:nAngs+p.nPsi,:)=[phi*ones(p.nPsi,1) theta*ones(p.nPsi,1) psis];
        nAngs=nAngs+p.nPsi;
    end;
end;

nAngs=size(angs,1);
p.nAngs=nAngs;
disp(p);

shiftVector=zeros(p.nPsi,2);
for i=1:p.nPsi
    shiftVector(i,:)=round(p.shiftMag*RotMatrix2((i-1)*2*pi/p.nPsi)*[1;0]); % shift along with psi
end;
shifts=repmat(shiftVector,nAngs/p.nPsi,1);
disp(p);

% -----------projections----------------
%     Making projections takes a long time, so we skip this if they have
%     already been computed.
makeProjections=~(exist('templates','var') && all(size(templates)==[p.imgSize p.imgSize nAngs]));

if makeProjections % do this if templates haven't already been calculated.
    fprintf(' making %d projections %d x %d\n',nAngs,p.imgSize,p.imgSize);
d1=struct;
d1.rlnAngleRot=angs(:,1);
d1.rlnAngleTilt=angs(:,2);
d1.rlnAnglePsi=angs(:,3);
d1.rlnOriginXAngst=shifts(:,1);
d1.rlnOriginYAngst=shifts(:,2);
templates=rlMakeRelionProjections(m,d1,[],p.pixA,p.dotCount);
%     rlAngs=angs;
%     rlAngs(:,3)=angs(:,3)-90; % fix a discrepancy
%     templates=rlMakeTemplates(rlAngs,m,dotCount);
disp(['Saving the templates at ' pwd '/templates.mat']);
save('templates.mat','templates','d1','p','-v7.3');
else
    disp('Using the existing templates.');
end; % if makeProjections

%%

nImgs=nAngs*nAmps; % we repeat the angles for each ampltude value.

nAmpMicrographs=ceil(nAngs/p.imgsPerMicrograph);
defStep=(p.defMax-p.defMin)/(nAmpMicrographs-1);
ctfs=zeros(p.imgSize,p.imgSize,nAmpMicrographs,'single');
lambda=EWavelength(p.kV);
Cs=2.7;
alpha=.02; % alpha for simulation

disp(' making the Optics struct');
    opt=struct;
opt.rlnOpticsGroupName={'opticsGroup1'};
    opt.rlnOpticsGroup=1;
    opt.rlnMicrographOriginalPixelSize=p.pixA;
    opt.rlnVoltage=p.kV;
    opt.rlnSphericalAberration=Cs;
    opt.rlnAmplitudeContrast=alpha;
    opt.rlnMicrographOriginalPixelSize=p.pixA;
    opt.rlnImagePixelSize=p.pixA;
    opt.rlnImageSize=p.imgSize; % we're setting the particle image size.
    opt.rlnImageDimensionality=2;

    disp(' making the CTFs');
defs=zeros(nAmpMicrographs,1);
for j=1:nAmpMicrographs
    defs(j)=p.defMin+(j-1)*defStep;
    ctfs(:,:,j)=CTF(p.imgSize,p.pixA,EWavelength(p.kV),defs(j),Cs,p.B,alpha);
end;

p.nImgs=nImgs;
p1.nImgs=nImgs;
disp(p1);

% To fill in
d.rlnImageName=cell(nImgs,1);
d.rlnMicrographName=cell(nImgs,1);
d.rlnDefocusU=zeros(nImgs,1);
d.rlnDefocusV=zeros(nImgs,1);
% Constant
d.rlnDefocusAngle=zeros(nImgs,1);
d.rlnCtfFigureOfMerit=ones(nImgs,1);
d.rlnOpticsGroup=ones(nImgs,1);
% Angles to fill in
d.rlnAngleRot=zeros(nImgs,1);
d.rlnAngleTilt=zeros(nImgs,1);
d.rlnAnglePsi=zeros(nImgs,1);
d.rlnOriginXAngst=zeros(nImgs,1);
d.rlnOriginYAngst=zeros(nImgs,1);

stack=zeros(p.imgSize,p.imgSize,nAngs*nAmps,'single');
%%
% Noise
disp(' making the noise...');

if p.useUniqueNoise
    nNoise=nImgs;
else
    nNoise=nAngs;
end;

if p.useWhiteNoise
    noise=single(randn(p.imgSize,p.imgSize,nNoise)); % white noise to add after CTF
else
    no.amp=.1;
    no.lFilt=.09;
    no.white=.04;
    p.noise=no;
    disp('Noise 1, before CTF');
    disp(no);
    
    no.amp=.05;
    no.lFilt=.083;
    no.white=.04;
    p.noise(2)=no;
    disp('Noise 2, after CTF');
    disp(no);
    
    nWhite=single(randn(p.imgSize,p.imgSize,nNoise,numel(p.noise))); % noise to be CTF-filtered
    
    noise=nWhite; % Get an array of the same size.
    
    for i=1:numel(p.noise)
        no=p.noise(i);
        noise(:,:,:,i)=no.amp*LorentzFilt(nWhite(:,:,:,i),no.lFilt,1)+no.white*nWhite(:,:,:,i);
    end;
end;

%%
maxRot=180/p.symmetry;
angMicrographIndex=floor((0:nAngs-1)'/p.imgsPerMicrograph)+1;
disp(' making the stack...');

%
for iAmp=1:nAmps
  for iAng=1:nAngs
    i=iAng+nAngs*(iAmp-1); % index of all images
    j1=angMicrographIndex(iAng);
    j=j1+(iAmp-1)*nAmpMicrographs; % index over all micrographs
    
    d.rlnMicrographName{i}=sprintf('%sm%04u%s',p.outDir,j,'.mrc');
    d.rlnImageName{i}=sprintf('%05u%s%s%s',i,'@',p.outDir,stackName);
    d.rlnDefocusU(i)=1e4*defs(j1);
    d.rlnDefocusV(i)=1e4*defs(j1);
%     rot=mod(angs(iAng,1)+maxRot,2*maxRot)-maxRot; % phi, restrict to +/- maxRot
    
    d.rlnAngleRot(i)=d1.rlnAngleRot(iAng);
    d.rlnAngleTilt(i)=d1.rlnAngleTilt(iAng);
    d.rlnAnglePsi(i)=d1.rlnAnglePsi(iAng);
    d.rlnOriginXAngst(i)=d1.rlnOriginXAngst(iAng);
    d.rlnOriginYAngst(i)=d1.rlnOriginYAngst(iAng);
    
    if p.useUniqueNoise
        ind=i; % different noise for each image
    else
        ind=iAng; % repeat the noise like we repeat the angles.
    end;
    %     Note that we are shifting after doing all the rotations.
    img=-p.amps(iAmp)*circshift(templates(:,:,iAng),-shifts(iAng,:));
    
    if p.useWhiteNoise
        stack(:,:,i)=real(ifftn(fftn(img).*ifftshift(ctfs(:,:,j1)))) ...
            +p.sigma*noise(:,:,ind);
    else
        stack(:,:,i)=real(ifftn(fftn(img+p.sigma*noise(:,:,ind,1)) ...
            .*ifftshift(ctfs(:,:,j1))))+p.sigma*noise(:,:,ind,2);
    end;
  end;
end;

%%
fullStackName=[p.outDir stackName];
disp(['Writing ' fullStackName]);
WriteMRC(stack,p.pixA,fullStackName);

%%
fullStarName=[p.outDir starName];
fullRefName=[p.outDir refName];
fullLogName=[p.outDir logName];

structNames={'data_optics'; 'data_particles'};
structs={opt; d};
disp(['Writing ' fullStarName]);
WriteStarFile(structNames,structs,fullStarName,'version 3');

disp(['Writing ' fullRefName]);
WriteMRC(m,p.pixA,fullRefName);

disp(['Writing ' fullLogName]);
WriteStructText(p,fullLogName);

disp('done.');

% figure(1);
% imagsar(stack);

