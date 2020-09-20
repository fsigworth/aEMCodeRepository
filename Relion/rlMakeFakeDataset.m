% rlMakeFakeDataset
% Comparing our angle assignments with Relion's
% rl_rot = -phi = -angs(i,1)
% rl_tilt=theta= angs(i,2)
% rl_psi = -psi-90 = -angs(i,3)-90;
% we apply shifts after rotating and projecting.
% The star file written out here corresponds to the data.star file produced in
% 3DRefine.

mapName='/Users/fred/aEMCodeRepository/AMPAR/KvMap.mat';
pa=fileparts(which('arGetRefVolumes'));
mapName=[pa '/KvMap.mat'];
outDir='/Users/fred/EMWork/Simulations/Relion/';
outDir='/home/siggpu/data/relion_sim/';
% stackName='SimStack03deg2.mrcs';
stackName='SimStack3amp15ang6.mrcs';
% starName='SimStack03deg2.star';
starName='SimStack3amp15ang6.star';
refName='Ref3Damp15.mrc';


amp=.015; % half amplitude
% % amp=.03 % particle signal
% amp=0  %% case for testing noise
sigma=11.925 % noise makes unity variance (empirical)
% sigma=0; %% case for noiseless projections
B=60       % ctf B factor

% makeProjections=~exist('templates','var');
makeProjections=true;

if makeProjections % do this if templates haven't already been calculated.
    
    s=load(mapName);  % gets s.map, s.pixA; map is 108^3 in size.
    % ShowSections(s.map);
    
    ds=1;
    imgSize=128;

    m=DownsampleGeneral(s.map,imgSize,1/ds);
    pixA=s.pixA*ds;
    
    symmetry=4;
    angStep=4;
%     angStep=6; %%%
    psiStep=90; % all angles are in degrees
    shiftMag=4;
    
    nTheta=round(180/angStep)+1
    dTheta=180/(nTheta-1);
    nPsi=ceil(360/psiStep)
    dPsi=360/nPsi;
    psis=(0:dPsi:360-dPsi)';
    maxNPhi=ceil(360/(angStep*symmetry))
    
    nAngs=0;
    angs=zeros(0,3);
    for i=1:nTheta
        theta=(i-1)*dTheta;
        nPhi=ceil(sind(theta)*maxNPhi); % Sample phi sparsely when sin(theta) is small.
        dPhi=360/(nPhi*symmetry);
        for j=1:nPhi
            phi=(j-1)*dPhi;
            angs(nAngs+1:nAngs+nPsi,:)=[phi*ones(nPsi,1) theta*ones(nPsi,1) psis];
            nAngs=nAngs+nPsi;
        end;
    end;
    
    nAngs=size(angs,1);
    nAngs
  
    shiftVector=zeros(nPsi,2);
    for i=1:nPsi
        shiftVector(i,:)=round(shiftMag*RotMatrix2((i-1)*2*pi/nPsi)*[1;0]); % shift along with psi
    end;
    shifts=repmat(shiftVector,nAngs/nPsi,1);
    
    templates=rlMakeTemplates(angs,m);
    
end; % if makeProjections
%%

imgsPerMicrograph=200;
kV=300;

nMicrographs=ceil(nAngs/imgsPerMicrograph)
defMin=1.5;  % Assign min, max defocus
defMax=3;
defStep=(defMax-defMin)/(nMicrographs-1);
ctfs=zeros(imgSize,imgSize,nMicrographs,'single');
lambda=EWavelength(kV);
Cs=2.7;
alpha=.02; % alpha for simulation
defs=zeros(nMicrographs,1);
for j=1:nMicrographs
    defs(j)=defMin+(j-1)*defStep;
    ctfs(:,:,j)=CTF(imgSize,pixA,EWavelength(300),defs(j),Cs,B,alpha);
end;

d=struct;
micrographIndex=floor((0:nAngs-1)'/imgsPerMicrograph)+1;
% To fill in
d.rlnImageName=cell(nAngs,1);
d.rlnMicrographName=cell(nAngs,1);
d.rlnDefocusU=zeros(nAngs,1);
d.rlnDefocusV=zeros(nAngs,1);
% Constant
d.rlnDefocusAngle=zeros(nAngs,1);
d.rlnVoltage=kV*ones(nAngs,1);
d.rlnAmplitudeContrast=.1*ones(nAngs,1); % set to usual Relion value.
d.rlnSphericalAberration=Cs*ones(nAngs,1);
d.rlnMagnification=5e4/pixA*ones(nAngs,1);
d.rlnDetectorPixelSize=5*ones(nAngs,1);
d.rlnCtfFigureOfMerit=ones(nAngs,1);
% Angles to fill in
d.rlnAngleRot=zeros(nAngs,1);
d.rlnAngleTilt=zeros(nAngs,1);
d.rlnAnglePsi=zeros(nAngs,1);
d.rlnOriginX=zeros(nAngs,1);
d.rlnOriginY=zeros(nAngs,1);

stack=zeros(imgSize,imgSize,nAngs,'single');

% Noise
n1=single(randn(imgSize,imgSize,nAngs)); % noise to be CTF-filtered
n2=single(randn(imgSize,imgSize,nAngs)); % noise added after CTF
noise1=.1*LorentzFilt(n1,.09,1)+.04*n1; % was .028
noise2=.05*LorentzFilt(n2,.083,1)+.04*n2; % was .026

maxRot=180/symmetry;

for i=1:nAngs
    j=micrographIndex(i);
    d.rlnMicrographName{i}=sprintf('m%04u%s',j,'.mrc');
    d.rlnImageName{i}=sprintf('%05u%s%s',i,'@',stackName);
    d.rlnDefocusU(i)=1e4*defs(j);
    d.rlnDefocusV(i)=1e4*defs(j);
    rot=mod(angs(i,1)+maxRot,2*maxRot)-maxRot; % phi, restrict to +/- maxRot
    d.rlnAngleRot(i)=rot;
    d.rlnAngleTilt(i)=angs(i,2);
    d.rlnAnglePsi(i)=angs(i,3);
    d.rlnOriginX(i)=shifts(i,1);
    d.rlnOriginY(i)=shifts(i,2);
    
%     Note that we are shifting after doing all the rotations.
    img=-amp*circshift(templates(:,:,i),-shifts(i,:))+sigma*noise1(:,:,i);
    stack(:,:,i)=real(ifftn(fftn(img) ...
        .*ifftshift(ctfs(:,:,j))))+sigma*noise2(:,:,i);
end;
%%

fullStackName=[outDir stackName];
fullStarName=[outDir starName];
fullRefName=[outDir refName];

disp(['Writing ' fullStackName]);
WriteMRC(stack,pixA,fullStackName);

disp(['Writing ' fullStarName]);
WriteStarFileStruct(d,'',fullStarName);

disp(['Writing ' fullRefName]);
WriteMRC(m,pixA,fullRefName);

disp('done.');
figure(1);
imagsar(stack);


return


%% Fit the background noise spectrum from a micrograph
nd=96;

micrName='/Users/fred/EMWork/Hideki/20181218/No5/sq04_1/Merged/Dec19_13.55.21ms.mrc';
[m,s]=ReadMRC(micrName);
imags(m);
mc=m(481:480+nd,756:755+nd);
figure(2);
subplot(221);
imags(mc);
drawnow;
spc=RadialPowerSpectrum(mc);
subplot(222);

semilogy(spc);

x=(1:48)';
f1=.01./(1+(x/9).^2)+.0008; % i.e. fc=9/96
f2=.0025./(1+(x/8).^2)+.0007; % i.e. fc=8/96
semilogy([spc f1 f2]);




