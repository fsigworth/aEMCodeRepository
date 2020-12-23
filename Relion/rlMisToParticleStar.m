% rlMisToParticleStar.m

infoDir='Info/';
starDir='';
starName='particles1.star';
useRawMicrograph=1; % Read the image file
useSubtractedImage=1;
setAllActive=1;
doPrint=1;

names=f2FindInfoFiles(infoDir);
ni=numel(names);
if ni<1
    disp(['No mi files found in ' infoDir]);
    return
end;

pts=struct;
j=0; % particle counter
imgSize=256; % just guessing
FlagRange=[16 32]; % flags for valid particles

nTotal=0;
ni=min(ni,50)
for i=1:ni
    miName=names{i};
    mi=ReadMiFile(miName);
    if i==1 % pick up optics parameters
        % make the optics structure
        opt=struct;
        opt.rlnOpticsGroupName={'opticsGroup1'};
        opt.rlnOpticsGroup=1;
        opt.rlnMicrographOriginalPixelSize=mi.pixA;
        opt.rlnVoltage=mi.kV;
        opt.rlnSphericalAberration=mi.ctf(1).Cs;
        opt.rlnAmplitudeContrast=mi.ctf(1).alpha;
        opt.rlnImagePixelSize=mi.pixA;
        opt.rlnImageSize=imgSize; % we're setting the particle image size.
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
        
        if useRawMicrograph
            if useSubtractedImage
                micName=[mi.procPath mi.baseFilename '_v.mrc'];
            else
                micName=[mi.imagePath imageFilenames{1}];
            end;
        else
            if useSubtractedImage
                micName=[mi.procPath mi.baseFilename 'mv.mrc'];
            else
                micName=[mi.procPath mi.baseFilename 'm.mrc'];
            end
        end;
        %     Accumulate the particles star
        istart=nTotal+1;
        iend=nTotal+nParts;
        pts.rlnMicrographName(istart:iend,1)={micName};
        pts.rlnCoordinateX(istart:iend,1)=xs;
        pts.rlnCoordinateY(istart:iend,1)=ys;
        nTotal=iend;
    end;
end;

% Fill in the constant fields
pts.rlnClassNumber(1:nTotal,1)=1;
pts.rlnAnglePsi(1:nTotal,1)=-999;
pts.rlnAutopickFigureOfMerit(1:nTotal,1)=-999;
pts.rlnOpticsGroup(1:nTotal,1)=1;

% Write the star file
partStarName=[starDir starName];
fStar=fopen(partStarName,'wt');
fprintf(fStar,'\n# version 30001\n');
WriteStarFileStruct(opt,'optics',fStar);
WriteStarFileStruct(pts,'particles',fStar);
fclose(fStar);

disp('Done.');

return
