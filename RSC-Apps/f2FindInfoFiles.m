function names=f2FindInfoFiles(infoDir)
% Searches the infoDir directory for *mi.txt files, and collects the names
% (infoDir path included) into the cell array names.  The infoDir argument
% is optional, in which case it is assumed to be Info/
if nargin<1
    infoDir='Info/';
else
    infoDir=AddSlash(infoDir);
end;
pattern='mi.txt';

d=dir(infoDir);

nNames=0;
names={};
for i=1:numel(d);
    if strndcmp(d(i).name,pattern)  % search the end of the name string
        nNames=nNames+1;
        names{nNames}=[infoDir d(i).name];
    end;
end;
disp([num2str(nNames) ' mi files found.']);
