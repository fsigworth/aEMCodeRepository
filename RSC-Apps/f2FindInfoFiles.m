function names=f2FindInfoFiles(infoDir,blanksOnly)
% Searches the infoDir directory for *mi.txt files, and collects the names
% (infoDir path included) into the cell array names.  The infoDir argument
% is optional, in which case it is assumed to be Info/
if nargin<1
    infoDir='Info/';
else
    infoDir=AddSlash(infoDir);
end;
if nargin<2
    blanksOnly=0;
end;
pattern='mi.txt';

d=dir(infoDir);

nNames=0;
names={};
for i=1:numel(d);
    if strndcmp(d(i).name,pattern)  % search the end of the name string
        theName=[infoDir d(i).name];
        if blanksOnly
            if mod(i,100)==0
                fprintf('.');
            end;
            mi=ReadMiFile(theName);
            nlog=numel(mi.log);
            if nlog>0 % skip if there's something in the log.
                disp([num2str([i nlog]) ' skip ' theName]);
                continue;
            end;
        end;
        nNames=nNames+1;
        names{nNames}=theName;
    end;
end;
disp([' ' num2str(nNames) ' mi files found.']);
