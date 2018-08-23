function infos=meScanHidekiFiles2(mi0, groupN, readIfPresent, rexp)
% function infos=meScanHidekiFiles2(mi0, groupN, readIfPresent, rexp)
% Look at the files in the given directory (concatenation of basePath and
% imagePath) and create a micrograph info structure mi for each set of exposure
% files.  The following fields are filled in:
% mi.imagePath
% mi.procPath
% mi.infoPath
% mi.baseFilename
% Save the structure in a file having the full name
% [serverdir infodir basename 'mi.mat']
% Processed images (e.g. merged images) will be written into infodir, e.g.
% [procdir basename 'm.mrc'].
% If the info file already exists, do not overwrite it but instead load its
% contents as mi.
% Then return a cell array of all the structures.

% serverdir='/Volumes/TetraData/EMWork/Hideki/110628/';
% imagedir='HS12slot4_GluR2-GFP/';
% procdir=[serverdir imagedir];
% infodir=imagedir;
% overwrite=0;

if nargin<4
rexp='\d+_\d+_\d+\.mrc';  % pick up files such as 001_002_91.mrc
rexp='\d+_.+_\d+\.mrc';  % pick up files such as 001_anytext00_91.mrc
rexp='.+_.+_\d+\.mrc';  % pick up files such as sq01_anytext00_91.mrc
end;

% Some other examples of regular expressions
% rexp='RawImage_\d*\.tif'  find RawImage_81.tif
% rexp='(?i)\d*(u|o)\d*\.(mrc|dm3)' ignore case, find e.g. 01u3000.MRC
% basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/120428/AMPA-R_lipo_with_carbon/';
% imagePath='micrograph/';
% procPath='Merged';
% infoPath='Info';
% overwrite=1;
% basePath
% 
% % if nargin < 1
% %     basePath='';
% % end;
% % if nargin<2
% %     imagePath='';
% % end;
% % if nargin<3
% %     procPath=imagePath;
% % end;
% % if nargin<4
% %     infoPath=procPath;
% % end;
% % if nargin<5
% %     overwrite=0;
% % end;
% 
% % Get a generic info structure
% mi0=meCreateMicrographInfoStruct11;
% 
% mi0.basePath=AddSlash(basePath);
% mi0.imagePath=AddSlash(imagePath);
% mi0.procPath=AddSlash(procPath);
% mi0.infoPath=AddSlash(infoPath);
searchDir=[mi0.basePath mi0.imagePath];
if ~DirectoryExists(searchDir)
    error(['Directory doesn''t exist: ' searchDir]);
elseif numel(dir(searchDir))<3
    error(['Directory is empty: ' searchDir]);
end;
names=FindFilenames([mi0.basePath mi0.imagePath],rexp);
nnm=numel(names);
disp([num2str(nnm) ' files found']);

infos={};
for i=1:numel(names)/groupN
    j=(i-1)*groupN+1;
    mi=mi0;  % copy the generic struct
    [pa mi.baseFilename ext]=fileparts(names{j});
    mi.imageFilenames{1}=names{j};
    for k=2:groupN
        mi.imageFilenames{k}=names{j-1+k};
    end;
    infoFilename=[mi.basePath mi.infoPath mi.baseFilename 'mi.mat'];
        if FileExists(infoFilename) && readIfPresent % the file exists already
            disp([infoFilename ' -loaded']);
            load(infoFilename);  % load the existing one.
        end;
    infos{i}=mi;
end;
    
