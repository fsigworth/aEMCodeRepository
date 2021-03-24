% rlMisToParticleStar.m
% Given a set of mi files and a micrographs_ctf.star file,
% create two particles.star files that can be used by
% Relion's particle extraction jobs. One is for unsubtracted and one is for
% subtracted particles. These files contain CTF parameters translated
% back from the mi files, and they contain two added fields, vesRadius and
% vesPsi.
%
% For unsubtracted particles we can use the raw micrograph or else a copy
% (possibly scaled) as Merged/*_u.mrc.
%
% Subtracted particles are either gotten from the padded, scaled micrograph
% Merged/*mv.mrc or from the unpadded micrograph Merged/*_v.mrc.
%
% We also create a micrographs_ctf_sub.star file that points to the subtracted
% micrographs, and if desired a micrographs_ctf_unsub.star file that points to our
% unsubtracted micrographs in the Merged/ folder.
%
% We can also write a vesicle star file which contains particle and vesicle
% coordinates, for predicting psi angles from geometry.



% ----Our picking data----
% First, use MiLoadAll to make an allMis.mat file containing all the mi file data.
% Then give the name here:
% allMisName='Picking_9/allMis9_intens+frac_7505.mat';
allMisName='Picking_9/allMis_holes_i2_ov_cls.mat';
%allMisName='allMis.mat';

% ----Inputs Micrograph star file
inMicStarName='CtfFind/job029/micrographs_ctf.star'; % Existing file to read

% ----Output star files
outStarDir='RSC9/';  % Place to put our particle star files
CheckAndMakeDir(outStarDir,1);

useMergedUnsubMicrograph=1; % Look for Merged/_u.mrc for unsub micrographs.
% We'll also write the new micrographs_ctf_unsub.star file.
usePaddedSubMicrograph=0; % Look for Merged/*mv.mrc for padded micrographs.
% Otherwise, look for Merged/*_v.mrc files.

outMicrographStarBasename='micrograph_ctf';
writeMicrographStarU=1;
writeMicrographStarV=1;

outParticleStarBasename='particles';
writeParticleStarU=1;
writeParticleStarV=1;

%   If desired, we write a subtracted file with vesicle coordinates and particle angles too.
outVesicleStarName=['ves_' outParticleStarBasename '_v.star'];
writeVesicleStar=0;
writeVesicleMat=0;

useGroupsFromMi=1; % Read the assigned group no. from mi.ok(20)
% OR ELSE just use an incrementing index of micrographs, with
minGroupParts=200; % minimun number of particles in a group

setParticlesActive=1; % ignore particle.picks(:,10) flag.
setMisActive=1; % ignore mi.active field.

checkForMicrographs=0;


if useMergedUnsubMicrograph
    unsubMicrographSuffix='_u.mrc';
end;

if usePaddedSubMicrograph
    disp('Not yet implemented: use padded sub micrograph. ');
    return
else
    subMicrographSuffix='_v.mrc'; % for image made in micrograph coordinates, instead of 'mv.mrc'
end;

disp(['Reading ' inMicStarName]);
[mcNames,mcDat]=ReadStarFile(inMicStarName);
% mcNames
opt=mcDat{1};
mic=mcDat{2};
nMics=numel(mic.rlnMicrographName);
disp([num2str(nMics) ' micrographs in star file.']);
disp(' ');
inMicPath=fileparts(mic.rlnMicrographName{1}); % pick up the path from the first one.
uMicNames=cell(nMics,1);
vMicNames=cell(nMics,1);
% %

disp(['Loading ' allMisName ' ...']);
load(allMisName); % Get allMis cell array
ni=numel(allMis);
disp([num2str(ni) ' mi files']);
disp(' ');
%%
pts=struct;
partSubMicName=cell(nMics,1);
partUnsubMicName=cell(nMics,1);

ves=struct; % structure for the vesicle info
outOpt=opt;
uMic=mic;
vMic=mic;       % copy the full micrograph star. We'll replace only the names

boxSize=256; % nominal starting size
FlagRange=[16 32]; % flags for valid particles
groupIndex=1; % if we're not reading from mi.ok
groupParts=0;
nTotal=0; % particle counter
j=0; % line counter
pSkip=0; % counter of micrographs with no particles
zSkip=0; % counter of micrographs with group=0
nSkip=0; % total micrographs skipped.
miSkip=0; % no. mi files skipped.
nBad=0;  % counter for mismatched micrograph names
namesMatched=1;

disp('Accumulating the structures. List: line; micrograph; particles; total particles.');
for i=1:ni
    mi=allMis{i};
    if ~isfield(mi,'opticsGroup')
        mi.opticsGroup=1;
    end;
    if i==1 % pick up optics parameters from the very first mi file, and
        %         put in a few more fields.
        nlOpt=numel(outOpt.rlnOpticsGroup);
        outOpt.rlnImagePixelSize=outOpt.rlnMicrographPixelSize; % copy the vector
        outOpt.rlnImageSize(1:nlOpt,1)=boxSize; % we're setting the default particle image size.
        outOpt.rlnImageDimensionality(1:nlOpt,1)=2;
    end;

            %         Get the micrograph names
        if useMergedUnsubMicrograph
            newUnsubMicName=[mi.procPath mi.baseFilename '_u.mrc'];
        else
            newUnsubMicName=[mi.imagePath mi.imageFilenames{1}];
        end;
        if usePaddedSubMicrograph
            disp('Not yet implemented, using padded sub micrographs.');
            return
        else
            newSubMicName=[mi.procPath mi.baseFilename '_v.mrc'];
        end;
        
        % We might have new files missing.
        
    [oldMicPath, oldMicBasename]=fileparts(mic.rlnMicrographName{i});
    nmLength=numel(oldMicBasename);
    [~,newVMicBasename]=fileparts(newSubMicName);
    [~,newUMicBasename]=fileparts(newUnsubMicName);
    if ~strncmp(oldMicBasename,newVMicBasename,nmLength) ...
          || ~strncmp(newUMicBasename,newVMicBasename,nmLength)% should match up to the end
      disp([ oldMicBasename '  ' newVMicBasename '  ' newUMicBasename]);
      return     
        miSkip=miSkip+1;
        continue; % we'll skip ahead
    end;
    
    if useGroupsFromMi
        groupIndex=mi.ok(20); % a zero groupIndex means a bad micrograph
        if groupIndex==0
            zSkip=zSkip+1;
            continue;
        end;
    end;
    if isfield(mi.particle,'picks') && numel(mi.particle.picks)>0 && groupIndex>0
        % ----- Accumulate the particle star data -----
        if size(mi.particle.picks,2)<10 || setParticlesActive % don't have the flag field
            flags=mi.particle.picks(:,3);
            mi.particle.picks(:,10)=(flags>=FlagRange(1)) & (flags <=FlagRange(2)); % all valid particles are active
        end;
        if setMisActive
            mi.active=true;
        end;
        active=(mi.particle.picks(:,10)>0) & mi.active;
        % ignore all particles when mi is not active.
        nParts=sum(active);
        
        if nParts<1
            pSkip=pSkip+1;
            nSkip=nSkip+1;
            continue;
        end;
        
        xs=mi.particle.picks(active,1);
        ys=mi.particle.picks(active,2);
        amps=mi.particle.picks(active,5);
        
        
        
        if mod(i,1000)==0
            disp(sprintf('%7d  %s %4d %8d',i,mi.baseFilename,nParts,nParts+nTotal));
        end;
        
        %     Accumulate the particles star
        istart=nTotal+1;
        iend=nTotal+nParts;
        pts.rlnCoordinateX(istart:iend,1)=xs;
        pts.rlnCoordinateY(istart:iend,1)=ys;
        pts.rlnAutopickFigureOfMerit(istart:iend,1)=amps;
        
        pts.rlnGroupName(istart:iend,1)={['group_' num2str(groupIndex)]};
        if ~useGroupsFromMi
            groupParts=groupParts+nParts;
            %  disp([groupParts groupIndex]);
            if groupParts>=minGroupParts
                groupLastParticle=iend;
                groupParts=0;
                groupIndex=groupIndex+1;
            end;
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
%         ves.vesMicrographName(istart:iend,1)={micName};
        ves.vesCenterX(istart:iend,1)=vxs;
        ves.vesCenterY(istart:iend,1)=vys;
        ves.vesR(istart:iend,1)=vrs;
        ves.vesPsi(istart:iend,1)=vpsis;
        ves.vesRsos(istart:iend,1)=rsos;
        ves.vesInds(istart:iend,1)=vInds;
        ves.ptlX(istart:iend,1)=xs;
        ves.ptlY(istart:iend,1)=ys;
        
        % Add the extra fields to the particle struct
        pts.vesicleRadius(istart:iend,1)=vrs;
        pts.vesiclePsi(istart:iend,1)=vpsis;
        
        if checkForMicrographs && nBad<20
            if ~exist(newUnsubMicName,2)
                disp([newUnsubMicName ' not found.']);
                nBad=nBad+1;
            end;
            if ~exist(newSubMicName,2)
                disp([newSubMicName ' not found.']);
                nBad=nBad+1;
            end;
        end;
        partUnsubMicName(istart:iend,1)={newUnsubMicName}; % The only fields that differ.
        partSubMicName(istart:iend,1)={newSubMicName};
        
        nTotal=iend;
    else
        nSkip=nSkip+1;
    end; % if particles

    vMicNames{i}=newSubMicName; % replace the micrograph name.
    uMicNames{i}=newUnsubMicName;
    if mic.rlnOpticsGroup(i)~=mi.opticsGroup % not one to one
        error(['Discrepancy in optics group indices at ' num2str(i)]);
    end;
end; % for loop over micrograph mi files

disp([num2str(nSkip) ' micrographs skipped, of ' num2str(ni)]);
disp(['Micrographs skipped with unassigned group: ' num2str(zSkip)]);
disp(['Micrographs skipped with group but no particles: ' num2str(pSkip)]);

if ~useGroupsFromMi
    % Make sure the last group has enough particles
    if groupParts<=minGroupParts && groupIndex>1
        groupNameCell=pts.rlnGroupName(groupLastParticle);
        pts.rlnGroupName(groupLastParticle+1:end)=groupNameCell;
    end;
end;

% --prepare the new micrograph.star structures
uMics=mic;
uMics.rlnMicrographName=uMicNames;
vMics=mic;
vMics.rlnMicrographName=vMicNames;

% --Prepare the particles.star structure
% Fill in the constant fields
pts.rlnClassNumber(1:nTotal,1)=1;
pts.rlnAnglePsi(1:nTotal,1)=-999;

uPts=pts;
uPts.rlnMicrographName=partUnsubMicName;
vPts=pts;
vPts.rlnMicrographName=partSubMicName;
ves.vesMicrographName=partSubMicName;


%%

if writeMicrographStarU
    % ----Write the sub micrographs star file----
    fullSubMicName=[outStarDir outMicrographStarBasename '_u.star'];
    disp(['Writing ' fullSubMicName]);
    fStar=fopen(fullSubMicName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(outOpt,'optics',fStar);
    WriteStarFileStruct(uMic,'micrographs',fStar);
    fclose(fStar);
end;

if writeMicrographStarV
    % ----Write the sub micrographs star file----
    fullSubMicName=[outStarDir outMicrographStarBasename '_v.star'];
    disp(['Writing ' fullSubMicName]);
    fStar=fopen(fullSubMicName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(outOpt,'optics',fStar);
    WriteStarFileStruct(vMic,'micrographs',fStar);
    fclose(fStar);
end;

% Write the particles star file
if writeParticleStarU
    outName=[outStarDir outParticleStarBasename '_u.star'];
    disp(['Writing ' outName '...']);
    fStar=fopen(outName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(outOpt,'optics',fStar);
    WriteStarFileStruct(vPts,'particles',fStar);
    fclose(fStar);
end;
%%
if writeParticleStarV
    outName=[outStarDir outParticleStarBasename '_v.star'];
    disp(['Writing ' outName '...']);
    fStar=fopen(outName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(outOpt,'optics',fStar);
    WriteStarFileStruct(uPts,'particles',fStar);
    fclose(fStar);
end;

%%
% Write the vesicle star file
if writeVesicleStar
    outName=[outStarDir outVesicleStarName];
    disp(['Writing ' outName '...']);
    fStar=fopen(outName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(outOpt,'optics',fStar);
    WriteStarFileStruct(ves,'vesicles',fStar);
    fclose(fStar);
end;

if writeVesicleMat
    [~,vnm]=fileparts(outVesicleStarName);
    outName=[outStarDir vnm '.mat'];
    disp(['Writing ' outName '...']);
    save(outName,'ves');
end;

disp('Done.');

