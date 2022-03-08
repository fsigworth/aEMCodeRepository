pwd% rlMisToParticleStar.m
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

boxSize=256;

% ----Our picking data----
% First, use MiLoadAll to make an allMis.mat file containing all the mi file data.
% Then give the name here:
% allMisName='Picking_9/allMis9_intens+frac_7505.mat';
%allMisName='Picking_9/allMis_holes_i2_ov_cls.mat';
allMisName='RSC2/allMis.mat';

% ----Inputs Micrograph star file
% inMicStarName='CtfFind/job003/micrographs_ctf.star'; % Existing file to read
inMicStarName='CtfFind/job007/micrographs_ctf.star'; % Existing file to  % 20211112 xchg dataset

% ----Output star files
outStarDir='RSC2/';  % Place to put our particle star files
CheckAndMakeDir(outStarDir,1);

useMergedUnsubMicrograph=0; % Look for Merged/_u.mrc for unsub micrographs.
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
outVesicleStarNameV=['ves_' outParticleStarBasename '_v.star'];
writeVesicleStar=1;
writeVesicleMat=1;

useGroupsFromMi=0; % Read the assigned group no. from mi.ok(20)
% OR ELSE just use an incrementing index of micrographs, with
minGroupParts=200; % minimun number of particles in a group
maxNMicrographs=inf; % limit the number of micrographs to consider

setParticlesActive=1; % ignore particle.picks(:,10) flag.
setMisActive=1; % ignore mi.active field.

checkForMicrographs=1;
skipLoadingMis=0;

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
%
nMics=numel(mic.rlnMicrographName);
if maxNMicrographs<nMics
    mic=TrimStructFields(mic,1,maxNMicrographs); %%%%%%%
    nMics=maxNMicrographs;
end;
disp([num2str(nMics) ' micrographs in star file.']);
disp(' ');

%inMicPath=fileparts(mic.rlnMicrographName{1}); % pick up the path from the first one.

uMicNames=cell(nMics,1);
vMicNames=cell(nMics,1);
% %

if ~skipLoadingMis
disp(['Loading ' allMisName ' ...']);
load(allMisName); % Get allMis cell array
end;
ni=numel(allMis);
disp([num2str(ni) ' mi files']);
disp(' ');
%
ni=nMics; %%%%%%%

pts=struct;
partSubMicName=cell(1,1);
partUnsubMicName=cell(1,1);

ves=struct; % structure for the vesicle info
outOpt=opt; % copy the optics from the micrograph.star
nOpt=numel(outOpt.rlnOpticsGroup);
    outOpt.rlnImageSize=boxSize*ones(nOpt,1);
    outOpt.rlnImageDimensionality=2*ones(nOpt,1);
    
boxSize=256; % nominal starting size
FlagRange=[16 32]; % flags for valid particles
groupIndex=1; % if we're not reading from mi.ok
groupParts=0;
nTotal=0; % particle counter
pSkip=0; % counter of micrographs with no particles
zSkip=0; % counter of micrographs with group=0
nSkip=0; % total micrographs skipped.
miSkip=0; % no. mi files skipped.

disp('Accumulating the structures. List: line; micrograph; particles; total particles.');
for i=1:ni
    mi=allMis{i};
%     if ~isfield(mi,'opticsGroup')
%         mi.opticsGroup=1;
%     end;
%     if i==1 % pick up optics parameters from the very first mi file, and
%         %   the input star file.
%         nlOpt=numel(outOpt.rlnOpticsGroup);
%         outOpt=opt;
% %         outOpt.rlnImagePixelSize=outOpt.rlnMicrographPixelSize; % copy the vector
% %         outOpt.rlnImageSize(1:nlOpt,1)=boxSize; % we're setting the default particle image size.
% %         outOpt.rlnImageDimensionality(1:nlOpt,1)=2;
% %         outOpt=rmfield(outOpt,'rlnMicrographPixelSize');
%     end;

            %    Get the micrograph name and find it in the inMicStar file.
            fullImageName=[mi.imagePath mi.imageFilenames{1}];
            match=strcmp(fullImageName,mic.rlnMicrographName);
            micStarLine=find(match);
            if numel(micStarLine)<1
                disp(['Image name couldn''t be matched. mi: ' fullImageName ' typical Star name: ' mic.rlnMicrographName{i}]);
                nSkip=nSkip+1;
                continue % skip this mi file.
            elseif numel(micStarLine)>1
                disp(['?duplicate micrograph names in star file? ' fullImageName]);
            end;
            micStarLine=micStarLine(1); % This is the index into the micrographs.star file.
            % Pick up parameters from the mic star file.
            opticsGroup=mic.rlnOpticsGroup(micStarLine);
            ctfMaxRes=mic.rlnCtfMaxResolution(micStarLine);
            ctfFOM=mic.rlnCtfFigureOfMerit(micStarLine);
            
            
        if useMergedUnsubMicrograph
            newUnsubMicName=[mi.procPath mi.baseFilename '_u.mrc'];
        else
            newUnsubMicName=[mi.imagePath mi.imageFilenames{1}]; % usual case.
        end;
        if usePaddedSubMicrograph
            disp('Not yet implemented, using padded sub micrographs.');
            return
        else
            newSubMicName=[mi.procPath mi.baseFilename '_v.mrc']; % usual case.
        end;
        
%         % We might have new files missing.
%         
%     [oldMicPath, oldMicBasename]=fileparts(mic.rlnMicrographName{j});
%     nmLength=numel(oldMicBasename);
%     [~,newVMicBasename]=fileparts(newSubMicName);
%     [~,newUMicBasename]=fileparts(newUnsubMicName);
%     if ~strncmp(oldMicBasename,newVMicBasename,nmLength) ...
%           || ~strncmp(newUMicBasename,newVMicBasename,nmLength)% should match up to the end
%       disp(['Discrepancy: ' oldMicBasename '  ' newVMicBasename '  ' newUMicBasename]);
%       return     
%     end;
    
    vMicNames{i}=newSubMicName; % replace the micrograph name.
    uMicNames{i}=newUnsubMicName;
    
    if useGroupsFromMi
        groupIndex=mi.ok(20); % a zero groupIndex means a bad micrograph
        if groupIndex==0
            nSkip=nSkip+1;
            continue;
        end;
    end;
        
    % Now starts code that is conditional on particle number.
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
            if checkForMicrographs
                if ~exist(newUnsubMicName,'file')
                    disp([newUnsubMicName ' not found.']); 
                end;
                if ~exist(newSubMicName,'file')
                    disp([newSubMicName ' not found.']);
                end;
            end;
            
        end;
        
        %     Accumulate the particles star
        istart=nTotal+1;
        iend=nTotal+nParts;
        pts.rlnCoordinateX(istart:iend,1)=xs;
        pts.rlnCoordinateY(istart:iend,1)=ys;
        pts.rlnClassNumber(istart:iend,1)=1;
        pts.rlnAutopickFigureOfMerit(istart:iend,1)=amps;
        pts.rlnAnglePsi(istart:iend,1)=-999; % We could assign these...
        pts.rlnOpticsGroup(istart:iend,1)=opticsGroup;
        % no image name, just micrograph name...
        partUnsubMicName(istart:iend,1)={newUnsubMicName}; % The only fields that differ.
        partSubMicName(istart:iend,1)={newSubMicName};       
%         
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
        %         pts.rlnAstigmatism(istart:iend,1)=-mi.ctf.deltadef*1e4;
 
        pts.rlnDefocusU(istart:iend,1)=(mi.ctf.defocus+mi.ctf.deltadef)*1e4;
        pts.rlnDefocusV(istart:iend,1)=(mi.ctf.defocus-mi.ctf.deltadef)*1e4;
        % not assigned: CtfBfactor, CtfMaxResolution, CtffigureOfMerit
        pts.rlnDefocusAngle(istart:iend,1)=mi.ctf.theta*180/pi;%         pts.rlnOpticsGroup(istart:iend,1)=mi.opticsGroup;
        pts.rlnCtfScalefactor(istart:iend,1)=1;
        pts.rlnPhaseShift(istart:iend,1)=0;
        pts.rlnCtfMaxResolution(istart:iend,1)=ctfMaxRes;
        pts.rlnCtfFigureOfMerit(istart:iend,1)=ctfFOM;
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
        
        % Alas, we can't add these extra fields to the particle struct
%        pts.vesicleRadius(istart:iend,1)=vrs;
%        pts.vesiclePsi(istart:iend,1)=vpsis;
        
        
        nTotal=iend;
    else
        nSkip=nSkip+1;
    end; % if particles

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
% pts.rlnClassNumber(1:nTotal,1)=1;
% pts.rlnAnglePsi(1:nTotal,1)=-999; % Alas!! this field causes Extract to hang.

uPts=pts;
uPts.rlnMicrographName=partUnsubMicName;
vPts=pts;
vPts.rlnMicrographName=partSubMicName;
ves.vesMicrographName=partSubMicName;


%%
% Write the micrograph star files
if writeMicrographStarU
    % ----Write the sub micrographs star file----
    fullSubMicName=[outStarDir outMicrographStarBasename '_u.star'];
    disp(['Writing ' fullSubMicName]);
    fStar=fopen(fullSubMicName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    % We just use the optics block from the input micrograph star file.
    WriteStarFileStruct(opt,'optics',fStar);
    WriteStarFileStruct(uMics,'micrographs',fStar);
    fclose(fStar);
end;

if writeMicrographStarV
    % ----Write the sub micrographs star file----
    fullSubMicName=[outStarDir outMicrographStarBasename '_v.star'];
    disp(['Writing ' fullSubMicName]);
    fStar=fopen(fullSubMicName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(opt,'optics',fStar);
    WriteStarFileStruct(vMics,'micrographs',fStar);
    fclose(fStar);
end;

% Write the particles star files
if writeParticleStarU
    outName=[outStarDir outParticleStarBasename '_u.star'];
    disp(['Writing ' outName]);
    fStar=fopen(outName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(outOpt,'optics',fStar);
    WriteStarFileStruct(uPts,'particles',fStar);
    fclose(fStar);
end;
%
if writeParticleStarV
    outName=[outStarDir outParticleStarBasename '_v.star'];
    disp(['Writing ' outName]);
    fStar=fopen(outName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(outOpt,'optics',fStar);
    WriteStarFileStruct(vPts,'particles',fStar);
    fclose(fStar);
end;

%
% Write the vesicle star files
if writeVesicleStar
    outName=[outStarDir outVesicleStarNameV];
    disp(['Writing ' outName]);
    fStar=fopen(outName,'wt');
    fprintf(fStar,'\n# version 30001\n');
    WriteStarFileStruct(outOpt,'optics',fStar);
    WriteStarFileStruct(ves,'vesicles',fStar);
    fclose(fStar);
end;

if writeVesicleMat
    [~,vnm]=fileparts(outVesicleStarNameV);
    outName=[outStarDir vnm '.mat'];
    disp(['Writing ' outName '...']);
    save(outName,'ves');
    disp(' ');
end;

disp('Done.');

