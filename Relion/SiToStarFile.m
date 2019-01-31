% SiToStarFile
% Create .star and .mrcs files from our si.mat and stack.mrc files.
batchMode=1;
siPath='Stack2/';
% in which case siPath and siName must be defined.

holeShots=7; % num,ber of micrographs to group together

dataName='images';
stackSuffixIn='stack.mrc';
stackSuffixOut='stack.mrcs';
stackUSuffixIn='ustack.mrc';
stackUSuffixOut='ustack.mrcs';
writeStacks=1;  % write out the modified stacks

modName=''; % extra string in name.
unsubToo=1;  % include the unsubtracted stack
cropSize=0;  % no cropping
setWeights=0; % change weights in mi file copy
doNormalize=1;  % run NormalizeImages
doGlobalNorm=0;
minSD=.1;  % criterion for including image at all
imgScale=1;
phasePlate=0;

shiftParticle=0; % shift the intracellular domain to center.
ctrOffsetA=24;  % IC domain center relative to old particle center
detPixel=5;     % pixel size of detector in um
ampContrast=0.1;  % if nonzero, override ctf(1).alpha value with this.

activeIndex=inf;
reverseContrast=0;

if ~batchMode
    % Put up a file selector for *si.mat,
disp('Getting an si file.');
[siName, siPath]=uigetfile('*si.mat','Select si file');
end;
if isnumeric(siPath)  % user has clicked Cancel
    return
end;
siPath=AddSlash(siPath);
disp(['Loading ' siPath siName]);
si=load([siPath siName]);
si=si.si;
nim=numel(si.miIndex);
afIndex=min(activeIndex,size(si.activeFlags,2));
activeFlags=si.activeFlags(:,afIndex);
fprintf('Active flag set: %d, %d active out of %d.\n',afIndex,sum(activeFlags),numel(activeFlags));

% Read the stack data
[pa,nm,ex]=fileparts(siName);
if numel(nm)<3
    error('si name is too short');
end;
inputStackName=[nm(1:end-2) stackSuffixIn];
inputStackUName=[nm(1:end-2) stackUSuffixIn];
disp(['Reading ' inputStackName]);
imgs0=imgScale*ReadMRC(inputStackName);

% Crop the stack
% if cropSize >0 && cropSize<size(imgs0,1)
%     disp(['Cropping to ' num2str(cropSize) ' pixels.']);
%     imgs0=Crop(imgs0,cropSize,1);
%     n=cropSize;
% else
    n=size(imgs0,1);
% end;
%% find normalization parameters: means, mesds

edgeMask=fuzzymask(n,2,n*0.45,0); % binary mask, dia=0.9*n
edgePix=edgeMask(:)>0;
nep=sum(edgePix);
imgEdgeVecs=zeros(nep,nim,'single');

for i=1:nim
    img=imgs0(:,:,i);
    imgEdgeVecs(:,i)=img(edgePix);
end;
means=mean(imgEdgeVecs,1)';
sds=std(imgEdgeVecs,1)';
% mem=mean(means);
% sdm=std(means);

%
af1=activeFlags;
sdWidth=2.5;
meWidth=3;

for ip=1:5
    if ip==1
        w=1.5;
    else
        w=1;
    end;
    mem=mean(means(af1));
    sdm=std(means(af1));
    mesds=mean(sds(af1));
    sdsds=std(sds(af1));
% disp(sum(af1));    
    af1=af1 & abs(sds-mesds)<sdWidth*w*sdsds & abs(means-mem)<meWidth*w*sdm;
end;

disp(['Overall std of images: ' num2str(mesds)]);
disp(['Original, activeFlags, final stack size: ' num2str([nim sum(activeFlags) sum(af1)])]);
figure(1);
subplot(221)
hist(means(af1),1000);
title('means (selected)');
subplot(222)
hist(means,1000);
title('all means');
subplot(223);
hist(sds(af1)/mesds,1000);
title('sds (norm, selected)');
subplot(224);
hist(sds,1000);
title('all sds');
%%


% Output file names
stackName=[nm(1:end-2) modName stackSuffixOut];
if unsubToo
    stackUName=[nm(1:end-2) modName stackUSuffixOut];
else
    stackUName='';
end;
starName=[nm(1:end-2) modName '.star'];
% disp(['Output files: ' starName '  ' stackName '  ' stackUName]);

%% Write the star file

disp(['Writing ' siPath starName]);
fi=fopen([siPath starName],'w');
% fi=1;

fprintf(fi,'%s %s\n#\n','# Original stack name: ',siName);
fprintf(fi,'data_%s\n',dataName);

fprintf(fi,'loop_\n');
fprintf(fi,'_rlnImageName\n');
fprintf(fi,'_rlnDefocusU\n');
fprintf(fi,'_rlnDefocusV\n');
fprintf(fi,'_rlnDefocusAngle\n');
fprintf(fi,'_rlnVoltage\n');
fprintf(fi,'_rlnAmplitudeContrast\n');
fprintf(fi,'_rlnSphericalAberration\n');
fprintf(fi,'_rlnMagnification\n');
fprintf(fi,'_rlnDetectorPixelSize\n');
fprintf(fi,'_rlnCtfFigureOfMerit\n');
% fprintf(fi,'_rlnGroupNumber\n');
fprintf(fi,'_rlnMicrographName\n');
if phasePlate
    fprintf(fi,'_rlnPhaseShift\n');
end;
if unsubToo
   fprintf(fi,'_rlnReconstructImageName\n');
end;
% fprintf(fi,'_rlnMicrographPart\n');

sinBetas=zeros(nim,1);
mag=detPixel*10000/si.pixA;
% ---------------need to fix unsub image stack selection and writing.---
for i=1:nim
%     q=normImgs(:,:,i);
%     s=std(q(:));
    if af1(i)  % ok particle, write a STAR line for it.
        miIndex=si.miIndex(i);
        mi=si.mi{miIndex};
        if any(setWeights)
            mi.weights=setWeights;
        end;
        ds=si.pixA/mi.pixA;
        if shiftParticle
            coords=mi.particle.picks(si.miParticle(i),:);
            rsoSign=sign(coords(7)-.5);  % +1 for rso, -1 for iso
            sinBeta=min(1,si.yClick(i)/(si.rVesicle(i)-rsoSign*si.mbnOffset));
            partShift=rsoSign*ctrOffsetA/si.pixA*sinBeta;
            yClick=si.yClick(i)-partShift;
%            normImgs(:,:,i)=Crop(circshift(normImgs(:,:,i),round([0 partShift])),cropSize);
            normImgs(:,:,i)=circshift(normImgs(:,:,i),round([0 partShift]));
        else
            alpha0=0;
        end;
        ctf=mi.ctf(1);
        defocus=ctf.defocus;
        deltadef=ctf.deltadef;        
        imgName=sprintf('%05d@%s',i,stackName);
        imgUName=sprintf('%05d@%s',i,stackUName);
        micrographPart=sprintf('%03d@%s',si.miParticle(i),mi.baseFilename);
        defU=round((defocus+deltadef)*1e4);
        defV=round((defocus-deltadef)*1e4);
        ang=mod(ctf.theta*180/pi-alpha0,180);
        kV=mi.kV;
        if ampContrast>0
            alpha=ampContrast;
        else
            alpha=ctf.alpha;
        end;
        if phasePlate
            if ~isfield(ctf,'phi')
                warning('No phi value in mi.ctf structure. No phase plate assumed.');
                ctf.phi=0;
                phasePlate=0;
            end;
            phaseShift=ctf.phi*180/pi;
        end;
        Cs=ctf.Cs;
        fom=si.sVesicle(i)*100;  % figure of merit is vesicle amplitude *100
        % Get the micrograph (group) name
        holeGroup=holeShots*floor((single(miIndex(i))-1)/holeShots)+1; % index of first image in group
%         if holeGroup>numel(si.mi)
%             error([num2str([si.miIndex(i) holeGroup]) ' holeGroup out of range.']);
%         end;
        gName=si.mi{holeGroup}.baseFilename;
            line=sprintf('%s %d %d %6.2f %g %5.3f %g %g %g %g %s',imgName,defU,defV,ang,kV,alpha,Cs,mag,detPixel,fom,gName);
        if phasePlate
            line=[line ' ' num2str(phaseShift)];
        end;
        if unsubToo
            line=[line ' ' imgUName];
        end;
%         line=[line ' ' micrographPart];
        fprintf(fi,'%s\n',line);        
    end;
end;
nim1=sum(af1);
fprintf(fi,'\n');
fclose(fi);
%% Normalize and write stacks
disp(['Normalization add,mul: ' num2str([mean(means) 1/mesds])]);
% normImgs=(imgs0-shiftdim(repmat(means,1,n,n),1))/mesds;
for i=1:nim
    normImgs(:,:,i)=(imgs0(:,:,i)-means(i))/mesds;
end;
if writeStacks
    disp(['Writing ' siPath stackName]);
    WriteMRC(normImgs,si.pixA,[siPath stackName]);
    
    if unsubToo
        disp(['Reading ' siPath inputStackUName]);
        imgs0=ReadMRC([siPath inputStackUName]);
        for i=1:nim
            normImgs(:,:,i)=(imgs0(:,:,i)-means(i))/mesds;
        end;
        disp(['Writing ' siPath stackUName]);
        WriteMRC(normImgs,si.pixA,[siPath stackUName]);
    end;
else
    disp('Not writing stacks.');
end;
%
%
% data_images
% loop_
% _rlnImageName
% _rlnDefocusU
% _rlnDefocusV
% _rlnDefocusAngle
% _rlnVoltage
% _rlnAmplitudeContrast
% _rlnSphericalAberration
% 000001@/lmb/home/scheres/data/VP7/all_images.mrcs 13538 13985 109.45 300 0.15 2
% 000002@/lmb/home/scheres/data/VP7/all_images.mrcs 13293 13796 109.45 300 0.15 2
% 000003@/lmb/home/scheres/data/VP7/all_images.mcrs 13626 14085 109.45 300 0.15 2