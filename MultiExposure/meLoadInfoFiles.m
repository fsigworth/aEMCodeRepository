function infos=meLoadInfoFiles(infoDir)
% function infos=meLoadInfoFiles(infodir)
% Look at the files in the given directory and load all the minfo files
% found there.
% Then return a cell array of all the structures.

if nargin < 1
    infoDir='Info/';
end;
infoDir=AddSlash(infoDir);

d=dir(infoDir);
nds=numel(d);
if nds<3
    error(['Bad directory: ' infoDir]);
end;
% We identify valid info file names as ones ending with
%   *mi.mat
infos={};
k=0;
% Scan the directory and pick up the name parts
for i=3:nds
    [pa base ext]=fileparts(d(i).name);
    len=numel(base);
    if len>4 && strcmp(ext,'.mat') && strcmp(base(len-1:len),'mi')
        load([infoDir base ext]);
        k=k+1;
        infos{k}=mi;
    end;
end;
