% SiToStarFile
% Create .star and .mrcs files from our si.mat and stack.mrc files.

dataName='images';
stackSuffix0='stack.mrc';
stackSuffix='stack.mrcs';
stackUSuffix0='ustack.mrc';
stackUSuffix='ustack.mrcs';
writeStacks=1;  % write out the modified stacks

modName='a1';
unsubToo=1;  % include the unsubtracted stack
cropSize=0;  % no cropping
setWeights=0; % change weights in mi file copy
doNormalize=1;  % run NormalizeImages
doGlobalNorm=0;
minSD=.5;  % criterion for including image
phasePlate=1;

shiftParticle=0; % shift the intracellular domain to center.
ctrOffsetA=24;  % IC domain center relative to old particle center
detPixel=5;     % pixel size in um
ampContrast=0;  % if nonzero, override ctf(1).alpha value with this.

activeIndex=inf;
reverseContrast=0;

% Put up a file selector for *si.mat,
disp('Getting an si file.');
[siName, siPath]=uigetfile('*si.mat','Select si file');
if ~isnumeric(siPath)  % user hasn't clicked Cancel
    cd(siPath);
end;
disp(siName);
si=load(siName);
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
inputStackName=[nm(1:end-2) stackSuffix0];
inputStackUName=[nm(1:end-2) stackUSuffix0];
disp(['Reading ' inputStackName]);
imgs0=ReadMRC(inputStackName);

% Crop the stack
% if cropSize >0 && cropSize<size(imgs0,1)
%     disp(['Cropping to ' num2str(cropSize) ' pixels.']);
%     imgs0=Crop(imgs0,cropSize,1);
%     n=cropSize;
% else
    n=size(imgs0,1);
% end;
%% purge extreme sds
af1=activeFlags;
sdWidth=3;
means=mean(reshape(imgs0,n*n,nim),1)'; % get whole-image means
for ip=1:3
    if ip==1
        sw=sdWidth*1.5;
    else
        sw=sdWidth;
    end;
    sds=std(reshape(imgs0,n*n,nim),1)';
    mesds=mean(sds(af1));
    sdsds=std(sds(af1));
    af1=activeFlags & sds>mesds-sw*sdsds & sds<mesds+sw*sdsds;
end;
disp(['Overall std of images: ' num2str(mesds)]);
disp(['Original, activeFlags, final stack size: ' num2str([nim sum(activeFlags) sum(af1)])]);
figure(1);
subplot(221)
hist(means(af1),1000);
title('means');
subplot(222)
hist(means,1000);
title('all means');
subplot(223);
hist(sds(af1),1000);
title('sds');
subplot(224);
hist(sds,1000);
title('all sds');

normImgs=(imgs0-shiftdim(repmat(means,1,n,n),1))/mesds;



% Output file names
stackName=[nm(1:end-2) modName stackSuffix];
if unsubToo
    stackUName=[nm(1:end-2) modName stackUSuffix];
else
    stackUName='';
end;
starName=[nm(1:end-2) modName '.star'];
disp(['Output files: ' starName ' ' stackName ' ' stackUName]);
%%

fi=fopen(starName,'w');
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
fprintf(fi,'_rlnGroupNumber\n');
if phasePlate
    fprintf(fi,'_rlnPhaseShift\n');
end;
if unsubToo
   fprintf(fi,'_rlnReconstructImageName\n');
end;
sinBetas=zeros(nim,1);
mag=detPixel*10000/si.pixA;
% ---------------need to fix unsub image stack selection and writing.---
for i=1:nim
    q=normImgs(:,:,i);
    s=std(q(:));
    if af1(i)  % ok particle
        mi=si.mi{si.miIndex(i)};
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
                error(['No phi value in mi.ctf structure']);
            end;
            phaseShift=ctf.phi*180/pi;
        end;
        Cs=ctf.Cs;
        fom=si.sVesicle(i)*100;  % figure of merit is vesicle amplitude *100
        gn=si.miIndex(i);
            line=sprintf('%s %d %d %6.2f %g %5.3f %g %g %g %g %g',imgName,defU,defV,ang,kV,alpha,Cs,mag,detPixel,fom,gn);
        if phasePlate
            line=[line ' ' num2str(phaseShift)];
        end;
        if unsubToo
            line=[line ' ' imgUName];
        end;
        fprintf(fi,'%s\n',line);        
    else
%        disp([num2str([i]) ' skipped.  SD= ' num2str(s)]);
    end;
end;
nim1=sum(af1);
fprintf(fi,'\n');
fclose(fi);
%%
% if doNormalize
%     [normImgs,vars,means]=NormalizeImages(imgs1,1,0,doGlobalNorm);  % global normalization
% %     if unsubToo
% %         normUImgs=imgsU0;
% %         for 
% %         normUIngs=(imgsU0-repmat(means,[imgSize(1:2) i]))./rep
%     disp(['Variance before normalization: ' num2str(median(vars))]);
% end;
% 
% if reverseContrast
%     normImgs=-normImgs;
% end;
% % eliminate outliers in std
% 
% upperThresh=3;  % beyond 3 sds is an outlier!
% lowerThresh=-3;
% normImgs(normImgs>upperThresh)=upperThresh;
% normImgs(normImgs<lowerThresh)=lowerThresh;

if writeStacks
    disp(['Writing ' stackName]);
    WriteMRC(normImgs,si.pixA,stackName);
    
    if unsubToo
        disp(['Reading ' inputStackUName]);
        normImgs=ReadMRC(inputStackUName)/mesds;
        disp(['Writing ' stackUName]);
        WriteMRC(normImgs,si.pixA,stackUName);
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