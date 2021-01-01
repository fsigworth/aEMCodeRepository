% rlMisToParticleStar.m
% Given a set of mi files, create a particles.star file that can be used by
% Relion's particle extraction job. This file contains CTF parameters translated
% back from the mi files. For unsubtracted particles we can use
% the raw micrograph. For subtracted micrographs we create our own
% micrographs_sub.star that points to our Merged directory and contains 
% just enough information to be used by Relion.
% At present we assume direct usage of the raw micrograph coordinates; that
% is, we don't use coordinates in padded micrographs.

infoDir='Info/';
starDir='';  % Place to put our star files
starName='particleSubC3groupMicNames.star';
micStarName='micrographs_sub.star';
useRawMicrograph=1; % Read the image file directly, use the raw coordinates.
useSubtractedImage=1;
setAllActive=1;
doPrint=1;
minGroupParts=100;

names=f2FindInfoFiles(infoDir);
ni=numel(names);
if ni<1
    disp(['No mi files found in ' infoDir]);
    return
end;

pts=struct;
mics=struct;
j=0; % particle counter
imgSize=256; % nominal starting size
FlagRange=[16 32]; % flags for valid particles
groupIndex=1;
groupParts=0;
nTotal=0;
%ni=min(ni,100)
for i=1:ni
    miName=names{i};
    mi=ReadMiFile(miName);
    if i==1 % pick up optics parameters from the very first mi file.
            % make the optics structure
        opt=struct;
        opt.rlnOpticsGroupName={'opticsGroup1'};
        opt.rlnOpticsGroup=1;
        opt.rlnMicrographOriginalPixelSize=mi.pixA;
        opt.rlnVoltage=mi.kV;
        opt.rlnSphericalAberration=mi.ctf(1).Cs;
        opt.rlnAmplitudeContrast=mi.ctf(1).alpha;
        micOpt=opt;         % optics paramtrs for the micrograph star file.
        micOpt.rlnMicrographPixelSize=mi.pixA;
        opt.rlnImagePixelSize=mi.pixA;
        opt.rlnImageSize=imgSize; % we're setting the default particle image size.
        opt.rlnImageDimensionality=2;
    end;
    
    if isfield(mi.particle,'picks') && numel(mi.particle.picks)>0
        if size(mi.particle.picks,2)<10 || setAllActive % don't have the flag field
            flags=mi.particle.picks(:,3);
            mi.particle.picks(:,10)=(flags>=FlagRange(1)) & (flags <=FlagRange(2)); % all valid particles are active
        end;
        active=mi.particle.picks(:,10)>0;
        nParts=sum(active);
        xs=mi.particle.picks(active,1);
        ys=mi.particle.picks(active,2);
        if doPrint && (i>1 && nParts>0)
            [~,nm,ex]=fileparts(miName);
            disp([nm ex '  ' num2str([nParts nParts+nTotal])]);
        end;
        
            if useSubtractedImage
                micName=[mi.procPath mi.baseFilename 'mv.mrc'];
            elseif useRawMicrograph
                micName=[mi.imagePath mi.imageFilenames{1}];
            else % use our processed, merged image.
                 micName=[mi.procPath mi.baseFilename 'mv.mrc']; 
            end
            % Note that 
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
        nTotal=iend;
    end;
    if useSubtractedImage % We make our own micrographs.star
        mics.rlnMicrographName{i}=[mi.procPath mi.baseFilename 'mv.mrc'];
        mics.rlnOpticsGroup(i)=1;
    end;
end;

% Make sure the last group is okay
    if groupParts<=minGroupParts && groupIndex>1
        groupNameCell=pts.rlnGroupName(groupLastParticle);
        pts.rlnGroupName(groupLastParticle+1:end)=groupNameCell;
    end;
% Fill in the constant fields
pts.rlnClassNumber(1:nTotal,1)=1;
pts.rlnAnglePsi(1:nTotal,1)=-999;
pts.rlnAutopickFigureOfMerit(1:nTotal,1)=-999;
pts.rlnOpticsGroup(1:nTotal,1)=1;

% Write the particles star file
partStarName=[starDir starName];
disp(['Writing ' partStarName]);
fStar=fopen(partStarName,'wt');
fprintf(fStar,'\n# version 30001\n');
WriteStarFileStruct(opt,'optics',fStar);
WriteStarFileStruct(pts,'particles',fStar);
fclose(fStar);
%%
if useSubtractedImage  % Write the micrographs star file
    micStarName=[starDir micStarName];
    disp(['Writing ' micStarName]);
fStar=fopen(micStarName,'wt');
fprintf(fStar,'\n# version 30001\n');
WriteStarFileStruct(micOpt,'optics',fStar);
WriteStarFileStruct(mics,'micrographs',fStar);
fclose(fStar);
end;

disp('Done.');

return
