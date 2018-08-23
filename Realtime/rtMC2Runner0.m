moviePath='movie_frames/';
movieExtension='.tif';
movieType='InTiff';
gainRefName='SuperRef.mrc';
imagePath='Micrograph/';
logPath='Log/';
pixA0=1.7/2;
frameDose=1;
kV=300;
ftBin=1;  % leave full size
gpus=0:3;
patches=[0 0];

MC2Exec='MotionCor2-1.0.0';
MC2Exec='/gpfs/ysm/apps/software/MotionCor2/1.0.5-gcccuda-2016.10/MotionCor2_Cuda8.0_v1.0.5';
MC2Exec='/gpfs/ysm/apps/software/MotionCor2/1.0.0-gcccuda-2016.10/MotionCor2-1.0.0';

doExec=1;

CheckAndMakeDir(imagePath);
CheckAndMakeDir(logPath);

names=GetFilenames(moviePath, movieExtension);
nf=numel(names);
disp([num2str(nf) ' files found.']);


for i=1:nf

disp([names{i} '  ']);
[pa mvBaseName ex]=fileparts(names{i});
spatches=num2str(patches);
sgpus=num2str(gpus);
sftBin=num2str(ftBin);

% Create the execution script
string=[MC2Exec ' -' movieType ' ' names{i}...
    ' -OutMRC ' [imagePath mvBaseName '.mrc'] ...
    ' -LogFile ' [logPath mvBaseName] '-' ...
    ' -Gain ' gainRefName ...
    ' -Patch ' spatches ' -Gpu ' sgpus ' -Kv ' num2str(kV) ...
    ' -FtBin ' sftBin ' -PixSize ' num2str(pixA0) ...
    ' -FmDose ' num2str(frameDose) ...
    ' >> ' [logPath 'Out.txt'] ];
disp(string);
if doExec
    system(string);
end;

end;

% %% Correct the micrograph scaling
% dwName=[imagePathPath mvBaseName '_DW.mrc']; % read the dose-weighted output
% [m,s]=ReadMRC(dwName);
% me=mean(m(:));
% m2=(m-me)*sqrt(nFrames);
% mOut=(Crop(m2,mi.imageSize)+me)*ds^2;  % pad the image, restore mean and scale up
% outName=[mi.imagePath mi.imageFilenames{1}];
% WriteMRC(mOut,s.pixA,outName);  % write it back out.
% disp(['Corrected micrograph written: ' outName]);
% 
% % delete the original files
% system(['rm ' mi.tempPath mvBaseName '*.mrc']);
% %%
% % Read the alignment shifts
% f=fopen([mi.tempPath mvBaseName '-0-Full.log']);
% header=fgetl(f);
% vals=cell2mat(textscan(f,'%f%f%f'));
% fclose(f);
% mi.frameShifts{1}=vals(:,2:3);
% 
% 
% 
% 
% plot(shifts);
% grid on;



function names=GetFilenames(path,pattern)
path=AddSlash(path);
d=dir(path);
nNames=0;
names={};
for i=3:numel(d)
    if any(strfind(d(i).name,pattern))
        names{end+1}=[path d(i).name];
    end;
end;
end