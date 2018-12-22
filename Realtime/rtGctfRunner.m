function [mi,epaVals,ctfImage,ctfVals]=rtGctfRunner(mi,mpars)
% function [mi,epaVals,ctfImage,ctfVals]=rtGctfRunner(mi,mpars);
%  We read the image file [mi.imagePath mi.imageFilename{1}] and determine
%  the ctf, based on the mi fields pixA and kV and the parameters given.  The
%  results are returned in the mi file, along with information for display.
% 
if nargin<2
    mpars=struct;
end;

% Parameter defaults:
pars.phasePlate=0;
pars.B=40;
pars.Cs=2.7;
pars.doPhaseFlip=0;
pars.minDefocus=0.2;
pars.maxDefocus=12;
pars.lowResolution=50;
pars.highResolution=4;
pars.alpha=.02;
pars.detectorPixelSize=5;
pars.gpuID=2;
pars.doDisplay=0;

pars=SetOptionValues(pars,mpars);
% mi.imageFilenames{1}=[mi.baseFilename 'ala.mrc'];
imageName=[mi.imagePath mi.imageFilenames{1}];
[~,baseImageName]=fileparts(mi.imageFilenames{1});

%GctfEx='Gctf.exe';
GctfEx='Gctf-v1.06_sm_30_cu8.0_x86_64';
if pars.phasePlate
    ppString=' --phase_shift_L 0 --phase_shift_H 180';
    pars.maxDefocus=1;
else
    ppString='';
end

% Construct the execution script and run it.
strings=cell(5,1);
strings{1}=[GctfEx ' --apix ' num2str(mi.pixA) ' --kV ' num2str(mi.kV) ' --cs ' num2str(pars.Cs) ' \'];
strings{2}=['--do_phase_flip ' num2str(pars.doPhaseFlip) ppString ' \'];
strings{3}=['--defL ' num2str(pars.minDefocus*1e4) ' --defH ' num2str(pars.maxDefocus*1e4) ...
    ' --resL ' num2str(pars.lowResolution) ' --resH ' num2str(pars.highResolution) ...
    ' --bfac ' num2str(pars.B) ' --ac ' num2str(pars.alpha) ' \'];
strings{4}=['--dstep ' num2str(pars.detectorPixelSize) ' --write_local_ctf 1' ' --do_EPA 1 \'];
strings{5}=[' --gid ' num2str(pars.gpuID) ' ' imageName ' >> ' mi.tempPath 'GctfOut.txt'];
%disp(strings{1});
exf=fopen('temp/ExecGctf.sh','w');
for i=1:numel(strings)
    disp(strings{i});
    fprintf(exf,'%s\n',strings{i});
end;
fclose(exf);
system('chmod a+x temp/ExecGctf.sh');
system('temp/ExecGctf.sh');

[epaVals,ctfVals,ctfImage]=rtReadGctfLogs(mi);

mi.ctf=struct;
mi.ctf.defocus=(ctfVals.Defocus_U + ctfVals.Defocus_V)/2e4;
mi.ctf.deltadef=(ctfVals.Defocus_U - ctfVals.Defocus_V)/2e4;
mi.ctf.theta=ctfVals.Angle*pi/180;
if pars.phasePlate
    mi.ctf.phi=ctfVals.Phase_shift*pi/180;
end;
mi.ctf.alpha=pars.alpha;
mi.ctf.lambda=EWavelength(mi.kV);
mi.ctf.Cs=pars.Cs;
mi.ctf.ampFactor=1;
mi.ctf.ccc=ctfVals.CCC;
mi.ctf.resLimit=ctfVals.RES_LIMIT;
mi.ctf.estB=ctfVals.B_FACTOR;
mi.ctf.B=pars.B;

if pars.doDisplay
    scl4=1/max(abs(epaVals.epaBkgSub));
    figure(1);
    subplot(2,2,2);
    imags(ctfImage);
    axis off;
    title(['Res limit: ' num2str(ctfVals.RES_LIMIT,3) 'Å']);
    subplot(2,1,2);    
    plot(1./epaVals.resolution,[epaVals.ctfSim.^2 scl4*epaVals.epaBkgSub epaVals.ccc 0*epaVals.ccc]);
    title(['\delta = ' num2str(mi.ctf.defocus,3) '\mum  \phi = ' num2str(ctfVals.Angle,3) '^o']);
end;

if pars.doPhaseFlip
    CheckAndMakeDir(mi.procPath);
    pfImageName=[mi.imagePath baseImageName '_pf.mrc'];
    pfOutName=[mi.procPath mi.baseFilename 'm.mrc'];
    pfOutNameSmall=[mi.procPath mi.baseFilename 'ms.mrc'];
    % The phase-flipped images is to be normalized.
    mpf=ReadMRC(pfImageName);
    me=mean(mpf(:));
    mpf2=(mpf/me-1);
    WriteMRC(mpf2,mi.pixA,pfOutName);
    mps=Downsample(mpf2,mi.imageSize/4);
    subplot(2,2,1);
    imags(mps);
    title(pfOutNameSmall,'interpreter','none');
    drawnow;
    WriteMRC(mps,mi.pixA*4,pfOutNameSmall);
%    system(['mv ' pfImageName ' ' pfOutName]);
end;

% Move all the log files from the image path to the temp path
if ~strcmp(mi.imagePath, mi.tempPath)
    str1=['mv ' mi.imagePath '*.log ' mi.tempPath];
    disp(str1);
    system(str1);
    % system(['mv ' mi.imagePath '*.ctf ' mi.tempPath]);
    % str2=['mv ' mi.tempPath 'micrographs_all_gctf.star ' mi.tempPath mi.baseFilename '.star'];
    % disp(str2);
    % system(str2);
end;
