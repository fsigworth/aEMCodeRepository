% rlMakeFakeMicrographs.m
% Derived from rlMakeFakeDataset
% Test of making picker outputs. The most convenient is to make a
% particles.star file for Relion extraction.
% 
% Meanwhile, comparing our angle assignments with Relion's,
% here is the result:
% rl_rot = -phi = -angs(i,1)
% rl_tilt=theta= angs(i,2)
% rl_psi = -psi-90 = -angs(i,3)-90;
% we apply shifts after rotating and projecting.

pa=fileparts(which('arGetRefVolumes'));
mapName=[pa '/KvMap.mat'];
% Output files
outDir=pwd; % Assume we're in the working directory
micDir='Micrographs2/';
starDir='Stars2/';
micBaseName='mic'

CheckAndMakeDir(micDir,1);
CheckAndMakeDir(starDir,1);

BFactor=60;
alpha=.05;
Cs=2.7;
kV=200;

nMicrographs=100
defRange=[3 6]
micSize=[4096 4096];
micBorder=80
minSpacing=100
ppm=500 % particlesPerMicrograph
maxParticlesPerMicrograph=prod(micSize-2*micBorder)/minSpacing^2

ds=1;
imgSize=128;
imgScale=.02; % 40 * 5e-4
symmetry=4;

shotSigma=3;
iceSigma=2*[1 .1]; % filtered+const
fcL=.05;
% iceSigma=11.925 % noise makes unity variance (empirical)


s=load(mapName);  % gets s.map, s.pixA; map is 108^3 in size.
% ShowSections(s.map);
map=DownsampleGeneral(s.map,imgSize,1/ds);
pixA=s.pixA*ds;

%     Make a linear ramp of defocus values
if nMicrographs>1
    ddef=diff(defRange)/(nMicrographs-1);
else
    ddef=0;
end;
defs=defRange(1)+ddef*(0:nMicrographs-1);
mics=zeros([micSize nMicrographs],'single');
pts=struct; % Fill up the particles star file.

rng(1);  % initialize random numbers
for i=1:nMicrographs
    psis=360*rand(ppm,1);
    thetas=acosd(1-2*rand(ppm,1));
    phis=360/symmetry*rand(ppm,1);
    templates=rlMakeTemplates([phis thetas psis],map);
    
    [xs,ys]=RandomTiling(ppm,micSize,minSpacing,micBorder);
    m=zeros(micSize,'single');
    for j=1:ppm
        m1=ExtractImage(templates(:,:,j),round([xs(j) ys(j)]),micSize,1);
        m=m+m1;
        % imags(m); drawnow;
    end;
    mics(:,:,i)=m; % store all the unfiltered micrographs
    imags(m);
    title(i);
    drawnow;
    
    % picks star file
    pk=struct;
    pk.rlnCoordinateX=xs;
    pk.rlnCoordinateY=ys;
    pk.rlnClassNumber=ones(ppm,1);
    pk.rlnAnglePsi=-999*ones(ppm,1);
    pk.rlnAutopickFigureOfMerit=-999*ones(ppm,1);
    starName=[micDir micBaseName num2str(i,'%03d') '_manualpick.star'];

    fpicks=fopen(starName,'wt');
    fprintf(fpicks,'\n# version 30001\n');
    WriteStarFileStruct(pk,'',fpicks);
    fclose(fpicks);
    
%     Write a .box file
    boxName=[micDir micBaseName num2str(i,'%03d') '_manualpick.box'];
    fbox=fopen(boxName,'wt');
    for j=1:ppm
        z=round([xs(j) ys(j)]-imgSize/2);
        fprintf(fbox,'%8u %8u %6u %6u\n',z(1),z(2),imgSize,imgSize);
    end;
    fclose(fbox);
    
    
    % accumulate for particles star file
    istart=(i-1)*ppm+1;
    iend=i*ppm;
    pts.rlnCoordinateX(istart:iend,1)=xs;
    pts.rlnCoordinateY(istart:iend,1)=ys;
    
end;

save([micDir 'mics.mat'],'mics','-v7.3');

%%
% make the optics structure
opt=struct;
for i=1 % I thought we had to make 2 lines to force a star table, but no.
    opt.rlnOpticsGroupName{i}=['opticsGroup' num2str(i)];
    opt.rlnOpticsGroup(i)=i;
    opt.rlnMicrographOriginalPixelSize(i)=pixA;
    opt.rlnVoltage(i)=kV;
    opt.rlnSphericalAberration(i)=Cs;
    opt.rlnAmplitudeContrast(i)=alpha;
    opt.rlnImagePixelSize(i)=pixA;
    opt.rlnImageSize(i)=imgSize; % we're setting the particle image size.
    opt.rlnImageDimensionality(i)=2;
end;


%Filter and write the micrographs
rng(2); % init the random number generator again
CheckAndMakeDir(micDir,1);
for i=1:nMicrographs
    
    c=CTF(micSize,pixA,EWavelength(kV),defs(i),Cs,BFactor,alpha);
    nIce=randn(micSize);
    iceNoise=iceSigma(1)*LorentzFilt(nIce,fcL/pixA)+iceSigma(2)*nIce;
    shotNoise=shotSigma*randn(micSize);
    m=mics(:,:,i);
    mc=real(ifftn(fftn(imgScale*m+iceNoise).*ifftshift(c)))+shotNoise;
    
    %     Write out the micrograph
    micName=[micDir micBaseName num2str(i,'%03d') '.mrc'];
    imags(GaussFilt(mc,.05/pixA));
    title(micName,'interpreter','none');
    drawnow;
    WriteMRC(mc,pixA,micName);
    %%
    
    %     Accumulate the particles star
    istart=(i-1)*ppm+1;
    iend=i*ppm;
    pts.rlnMicrographName(istart:iend,1)={micName};
    
end;
%%
nim=numel(pts.rlnCoordinateX);

% % Create the fake stack
% stackName='fakeStack.mrcs';
% stk=zeros(16,16,nim,'single');
%    for k=1:nim
%         pts.rlnImageName{k,1}=[num2str(k,'%04d') '@' stackName];
%     end;
% WriteMRC(stk,pixA,stackName);

% -----Make the particles.star file for Relion extraction
% We've set up the optics struct.
% And we have already set these fields:
% pts.rlnMicrographName
% pts.rlnCoordinateX
% pts.rlnCoordinateY
% now for the rest:
pts.rlnClassNumber(1:nim,1)=1;
pts.rlnAnglePsi(1:nim,1)=-999;
pts.rlnAutopickFigureOfMerit(1:nim,1)=-999;
pts.rlnOpticsGroup(1:nim,1)=1;

partStarName=[starDir 'particles.star'];
fStar=fopen(partStarName,'wt');
fprintf(fStar,'\n# version 30001\n');
WriteStarFileStruct(opt,'optics',fStar);
WriteStarFileStruct(pts,'particles',fStar);
fclose(fStar);

WriteMRC(map,pixA,'KvRef.mrc');

disp('Done.');

return

