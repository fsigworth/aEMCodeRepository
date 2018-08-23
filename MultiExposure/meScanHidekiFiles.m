function infos=meScanHidekiFiles(serverdir,imagedir, procdir, infodir, overwrite, seqMode)
% function infos=meScanLeginonFiles(serverdir,localdir,infodir, overwrite, seqMode)
% Look at the files in the given directory (concatenation of serverdir and
% imagedir) and create a micrograph info structure mi for each set of exposure
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


if nargin < 1
    serverdir='';
end;
if nargin<2
    localdir='';
end;
if nargin<3
    procdir=[serverdir localdir];
end;
if nargin<4
    infodir=procdir;
end;
if nargin<5
    overwrite=0;
end;

% Make sure serverdir ends with a slash.
len=numel(serverdir);
if len>0 && serverdir(len)~='/'
    serverdir=[serverdir '/'];
end;

% scan the current directory, looking for Hideki's files.  For each
% set of exposures, create a micrograph info structure, store it as an
% XXXmi.mat file, and also add it to the infos array of structures.

d=dir([serverdir imagedir]);
ind=0;
infos={};

for i=1:numel(d)
    [names basename defoci ind]=FindHidekiFileSet(d,ind+1,seqMode);
    if ind<=0
        break
    end;
    if numel(basename)>0
        infoFilename=[serverdir infodir basename 'mi.mat'];
        if FileExists(infoFilename) && ~overwrite % the file exists already
            disp([infoFilename ' -exists']);
            load(infoFilename);  % load the existing one.
            mi=meUpdateInfoStruct(mi);  % Convert to latest version
        else    % construct the mi struct from scratch
            mi=meCreateMicrographInfoStruct10;
            mi.baseFilename=basename;
            mi.imagePath=imagedir;
            mi.procPath=procdir;
            mi.infoPath=infodir;
            
            % Find out how many exposures we have, including extension '.mrc'
            mi.imageFilenames=names;
            for j=1:numel(defoci)
                mi.ctf(j).defocus=defoci(j);
            end;
            save(infoFilename,'mi');
            disp([infoFilename ' -created'])
        end;
        infos{i}=mi;
    end;
end;
