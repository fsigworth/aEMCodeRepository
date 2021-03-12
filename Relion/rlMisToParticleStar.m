% rlMisToParticleStar.m
% Given a set of mi files, create a particles.star file that can be used by
% Relion's particle extraction job. This file contains CTF parameters translated
% back from the mi files. For unsubtracted particles we can use
% the raw micrograph. For subtracted micrographs we create our own
% micrographs_sub.star that points to our Merged directory and contains
% just enough information to be used by Relion.
% At present we assume direct usage of the raw micrograph coordinates; that
% is, we don't use coordinates in padded micrographs.


% ----Our picking data----
% First, use MiLoadAll to make an allMis.mat file containing all the mi file data.
% Then give the name here:
allMisName='allMis.mat';

% ----Micrograph star files
micStarName='CtfFind/job003/micrographs_ctf.star';
useRawMicrograph=1; % Read unpadded unsub images rather than from the Merged directory
useSubtractedMicrograph=1; % Use the subtracted micrograph name in the particles file.
    % (The subtracted micrograph is assumed to be in the Merged/ folder.)
writeSubMicrographsStar=1; % write a new star file pointing to the sub micrographs?
subMicStarName='CtfFind/job003/micrographs_sub_ctf.star'; % New star file to write

% -----Particle and Vesicle info files to write-----
outStarDir='RSC/';  % Place to put our particle star files
    CheckAndMakeDir(outStarDir,1);
outParticleStarName='particles2_sub.star';
outVesicleStarName=['ves_' outParticleStarName];
writeParticleStar=1;
writeVesicleStar=1;
writeVesicleMat=1; % Instead of writing a long .star file, save as a Matlab .mat


setParticlesActive=1; % ignore particle.picks(:,10) flag.
doPrint=1;
minGroupParts=200; % minimun number of particles in a group


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
disp([num2str(numel(mic.rlnMicrographName)) ' micrographs in star file.']);
% %

disp(['Loading ' allMisName ' ...']);
load(allMisName); % Get allMis cell array
ni=numel(allMis);
disp([num2str(ni) ' mi files']);

%%

pts=struct;
ves=struct; % structure for the vesicle info
sOpt=opt;      % copy the optics info to the sub micrograph structure
sMic=mic;       % copy the full micrograph star. We'll replace only the names

boxSize=256; % nominal starting size
FlagRange=[16 32]; % flags for valid particles
groupIndex=1;
groupParts=0;
nTotal=0; % particle counter
j=0; % line counter

disp('Accumulating the structures. List: line; micrograph; particles; total particles.');
for i=1:ni
    %     miName=names{i};
    %     mi=ReadMiFile(miName);
    mi=allMis{i};
    if i==1 % pick up optics parameters from the very first mi file, and
        %         put in a few more fields.
        nlOpt=numel(opt.rlnOpticsGroup);
        opt.rlnImagePixelSize=opt.rlnMicrographPixelSize; % copy the vector
        opt.rlnImageSize(1:nlOpt,1)=boxSize; % we're setting the default particle image size.
        opt.rlnImageDimensionality(1:nlOpt,1)=2;
    end;
    if isfield(mi.particle,'picks') && numel(mi.particle.picks)>0
        % ----- Accumulate the particle star data -----
        if size(mi.particle.picks,2)<10 || setParticlesActive % don't have the flag field
            flags=mi.particle.picks(:,3);
            mi.particle.picks(:,10)=(flags>=FlagRange(1)) & (flags <=FlagRange(2)); % all valid particles are active
        end;
        if ~isfield(mi,'active')
            mi.active=1;
        end;
        if ~isfield(mi,'opticsGroup')
            mi.opticsGroup=1;
        end;
        active=(mi.particle.picks(:,10)>0) & mi.active; % ignore all particles when mi is not active.
        nParts=sum(active);
        
        if nParts<1
            continue;
        end;
        
        xs=mi.particle.picks(active,1);
        ys=mi.particle.picks(active,2);
        amps=mi.particle.picks(active,5);

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
        pts.rlnAutopickFigureOfMerit(istart:iend,1)=amps;

        pts.rlnGroupName(istart:iend,1)={['group_' num2str(groupIndex)]};
        groupParts=groupParts+nParts;
        %  disp([groupParts groupIndex]);
        if groupParts>=minGroupParts
            groupLastParticle=iend;
            groupParts=0;
            groupIndex=groupIndex+1;
        end;
        
        % For reference, this is how we got the mi.ctf parameters from the original star files:
        % mi.ctf.defocus=(mic.rlnDefocusU(iLine)+mic.rlnDefocusV(iLine))/2e4;
        % mi.ctf.deltadef=(mic.rlnDefocusU(iLine)-mic.rlnDefocusV(iLine))/2e4;
        % mi.ctf.theta=mic.rlnDefocusAngle(iLine)*pi/180;
        
        pts.rlnDefocusU(istart:iend,1)=(mi.ctf.defocus+mi.ctf.deltadef)*1e4;
        pts.rlnDefocusV(istart:iend,1)=(mi.ctf.defocus-mi.ctf.deltadef)*1e4;
        %         pts.rlnAstigmatism(istart:iend,1)=-mi.ctf.deltadef*1e4;
        pts.rlnDefocusAngle(istart:iend,1)=mi.ctf.theta*180/pi;
        pts.rlnOpticsGroup(istart:iend,1)=mi.opticsGroup;
                
        % ----- Accumulate the vesicle star -----
        rsos=mi.particle.picks(active,7); % rso flags
        vInds=mi.particle.picks(active,4); %vesicle indices
            % handle particles with no vesicle index
        vesOk=vInds>0;
        if any(~vesOk)
            disp(['Bad vesicle in image ' num2str(i) '  ' mi.baseFilename])
        end;
        vIndsOk=vInds(vesOk);
        vxs=zeros(nParts,1,'single');
        vxs(vesOk)=mi.vesicle.x(vIndsOk);
        vys=zeros(nParts,1,'single');
        vys(vesOk)=mi.vesicle.y(vIndsOk);
        vrs=zeros(nParts,1,'single');
        vrs(vesOk)=real(mi.vesicle.r(vIndsOk,1));
        vpsis=atan2d(ys-vys,xs-vxs);
        
        ves.vesMicrographName(istart:iend,1)={micName};
        ves.vesCenterX(istart:iend,1)=vxs;
        ves.vesCenterY(istart:iend,1)=vys;
        ves.vesR(istart:iend,1)=vrs;
        ves.vesPsi(istart:iend,1)=vpsis;
        ves.vesRsos(istart:iend,1)=rsos;
        ves.vesInds(istart:iend,1)=vInds;
        ves.ptlX(istart:iend,1)=xs;
        ves.ptlY(istart:iend,1)=ys;
        
                nTotal=iend;            
    else
        continue;
    end; % if particles

    if useSubtractedMicrograph % We make our own micrographs.star
        %             We're assuming here a one-to-one correspondence between mis
        %             and lines of the micrograph_ctf file.
        j=i; % line index. Note that we don't overwrite the original mic names where there
             %  are no particles!
        sMic.rlnMicrographName{j}=subMicName;
        if mic.rlnOpticsGroup(j)~=mi.opticsGroup % not one to one
            error(['Discrepancy in micrograph indices at ' num2str(j)]);
        end;
    end;
end; % for loop over micrograph mi files

% Make sure the last group is okay
if groupParts<=minGroupParts && groupIndex>1
    groupNameCell=pts.rlnGroupName(groupLastParticle);
    pts.rlnGroupName(groupLastParticle+1:end)=groupNameCell;
end;

% --Prepare the particles.star structure
% Fill in the constant fields
pts.rlnClassNumber(1:nTotal,1)=1;
pts.rlnAnglePsi(1:nTotal,1)=-999;

% Write the particles star file
if writeParticleStar
    outName=[outStarDir outParticleStarName];
    disp(['Writing ' outName '...']);
    fStar=fopen(outName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(opt,'optics',fStar);
    WriteStarFileStruct(pts,'particles',fStar);
    fclose(fStar);
end;
%
% Write the vesicle star file
if writeVesicleStar
    outName=[outStarDir outVesicleStarName];
    disp(['Writing ' outName '...']);
    fStar=fopen(outName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(opt,'optics',fStar);
    WriteStarFileStruct(pts,'vesicles',fStar);
    fclose(fStar);
end;

if writeVesicleMat
    [~,vnm]=fileparts(outVesicleStarName);
    outName=[outStarDir vnm '.mat'];
    disp(['Writing ' outName '...']);
    save(outName,'ves');
end;
%
if useSubtractedMicrograph && writeSubMicrographsStar
    % ----Write the sub micrographs star file----
    disp(['Writing ' subMicStarName '...']);
    fStar=fopen(subMicStarName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(sOpt,'optics',fStar);
    WriteStarFileStruct(sMic,'micrographs',fStar);
    fclose(fStar);
end;

disp('Done.');

