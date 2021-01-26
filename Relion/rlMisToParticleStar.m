% rlMisToParticleStar.m
% Given a set of mi files, create a particles.star file that can be used by
% Relion's particle extraction job. This file contains CTF parameters translated
% back from the mi files. For unsubtracted particles we can use
% the raw micrograph. For subtracted micrographs we create our own
% micrographs_sub.star that points to our Merged directory and contains
% just enough information to be used by Relion.
% At present we assume direct usage of the raw micrograph coordinates; that
% is, we don't use coordinates in padded micrographs.
allMisName='allMis8_sel.mat';
micStarName='CtfFind/job026/micrographs_ctf.star';
subMicStarName='CtfFind/job026/micrographs_sub_ctf.star';
% infoDir='Info/';
outStarDir='';  % Place to put our star files
outStarName='particleAllSubMicNames.star';
outStarName='particleAllMicNames.star';

useRawMicrograph=1; % Read unpadded images
useSubtractedMicrograph=0; % Use the subtracted micrograph name in the particles file.
% Also create the subtracted micrograph star file.

setParticlesActive=1; % ignore particle.picks(:,10) flag.
doPrint=1;
minGroupParts=100; % minimun number of particles in a group


if useRawMicrograph
    subMicrographSuffix='_v.mrc'; % for image made in micrograph coordinates, instead of 'mv.mrc'
else
    subMicrographSuffix='mv.mrc';
end;

% names=f2FindInfoFiles(infoDir);
% ni=numel(names);
% if ni<1
%     disp(['No mi files found in ' infoDir]);
%     return
% end;
%

disp(['Reading ' micStarName]);
[mcNames,mcDat]=ReadStarFile(micStarName);
mcNames
opt=mcDat{1};
mic=mcDat{2};
sOpt=opt;      % copy the optics info
sMic=mic;       % copy the subtracted micrograph star structure
disp([num2str(numel(mic.rlnMicrographName)) ' micrographs in star file.']);
% %

disp(['Loading ' allMisName ' ...']);
load(allMisName); % Get allMis cell array
disp(' done.');
ni=numel(allMis);
disp([num2str(ni) ' mi files']);

%%
pts=struct;
% mics=struct;
j=0; % particle counter
imgSize=256; % nominal starting size
FlagRange=[16 32]; % flags for valid particles
groupIndex=1;
groupParts=0;
nTotal=0;
%ni=min(ni,100)
for i=1:ni
    %     miName=names{i};
    %     mi=ReadMiFile(miName);
    mi=allMis{i};
    if i==1 % pick up optics parameters from the very first mi file, and
        %         put in a few more fields.
        nlOpt=numel(opt.rlnOpticsGroup);
        opt.rlnImagePixelSize=opt.rlnMicrographPixelSize;
        opt.rlnImageSize(1:nlOpt,1)=imgSize; % we're setting the default particle image size.
        opt.rlnImageDimensionality(1:nlOpt,1)=2;
    end;
    if isfield(mi.particle,'picks') && numel(mi.particle.picks)>0
        if size(mi.particle.picks,2)<10 || setParticlesActive % don't have the flag field
            flags=mi.particle.picks(:,3);
            mi.particle.picks(:,10)=(flags>=FlagRange(1)) & (flags <=FlagRange(2)); % all valid particles are active
        end;
        active=(mi.particle.picks(:,10)>0) & mi.active; % ignore all particles when mi is not active.
        nParts=sum(active);
        xs=mi.particle.picks(active,1);
        ys=mi.particle.picks(active,2);
        subMicName=[mi.procPath mi.baseFilename subMicrographSuffix];
        if useSubtractedMicrograph
            micName=subMicName;
        else
            if useRawMicrograph
                micName=[mi.imagePath mi.imageFilenames{1}];
            else
                micName=[mi.procPath mi.baseFilename 'm.mrc'];
            end;
        end;

        if doPrint && mod(i,1000)==0
            disp(sprintf('%7d  %s %4d %8d',i,mi.baseFilename,nParts,nParts+nTotal));
        end;
        
        %     Accumulate the particles star
        istart=nTotal+1;
        iend=nTotal+nParts;
        pts.rlnMicrographName(istart:iend,1)={micName};
        pts.rlnCoordinateX(istart:iend,1)=xs;
        pts.rlnCoordinateY(istart:iend,1)=ys;
        pts.rlnGroupName(istart:iend,1)={['group_' num2str(groupIndex)]};
        groupParts=groupParts+nParts;
        %  disp([groupParts groupIndex]);
        if groupParts>=minGroupParts
            groupLastParticle=iend;
            groupParts=0;
            groupIndex=groupIndex+1;
        end;
        
        % This is how we got the mi.ctf parameters from the original star files:
        % mi.ctf.defocus=(mic.rlnDefocusU(iLine)+mic.rlnDefocusV(iLine))/2e4;
        % mi.ctf.deltadef=(mic.rlnDefocusU(iLine)-mic.rlnDefocusV(iLine))/2e4;
        % mi.ctf.theta=mic.rlnDefocusAngle(iLine)*pi/180;
        
        pts.rlnDefocusU(istart:iend,1)=(mi.ctf.defocus+mi.ctf.deltadef)*1e4;
        pts.rlnDefocusV(istart:iend,1)=(mi.ctf.defocus-mi.ctf.deltadef)*1e4;
        %         pts.rlnAstigmatism(istart:iend,1)=-mi.ctf.deltadef*1e4;
        pts.rlnDefocusAngle(istart:iend,1)=mi.ctf.theta*180/pi;
        pts.rlnOpticsGroup(istart:iend,1)=mi.opticsGroup;
        nTotal=iend;
    end; % if particles
    if useSubtractedMicrograph % We make our own micrographs.star
        %             We're assuming here a one-to-one correspondence between mis
        %             and lines of the micrograph_ctf file.
        sMic.rlnMicrographName{i}=subMicName;
        if mic.rlnOpticsGroup(i)~=mi.opticsGroup % not one to one
            error(['Discrepancy in micrograph indices at ' num2str(i)]);
        end;
    end;
end; % for

% Make sure the last group is okay
if groupParts<=minGroupParts && groupIndex>1
    groupNameCell=pts.rlnGroupName(groupLastParticle);
    pts.rlnGroupName(groupLastParticle+1:end)=groupNameCell;
end;

% --Prepare the particles.star structure
% Fill in the constant fields
pts.rlnClassNumber(1:nTotal,1)=1;
pts.rlnAnglePsi(1:nTotal,1)=-999;
pts.rlnAutopickFigureOfMerit(1:nTotal,1)=-999;

% Write the particles star file
partStarName=[outStarDir outStarName];
disp(['Writing ' partStarName]);
fStar=fopen(partStarName,'wt');
fprintf(fStar,'\n# version 30001\n');
WriteStarFileStruct(opt,'optics',fStar);
WriteStarFileStruct(pts,'particles',fStar);
fclose(fStar);
%
%%
if useSubtractedMicrograph
    % ----Write the sub micrographs star file----
    fullSubMicName=[outStarDir subMicStarName];
    disp(['Writing ' fullSubMicName]);
    fStar=fopen(fullSubMicName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(sOpt,'optics',fStar);
    WriteStarFileStruct(sMic,'micrographs',fStar);
    fclose(fStar);
end;

disp('Done.');

