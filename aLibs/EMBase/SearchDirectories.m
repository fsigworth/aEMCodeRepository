function fNames=SearchDirectories(inDir,fNames,pars)
% Recursive search for files having the correct extensions or endings
% or full names. Builds up the
% cell array of complete filenames fNames. Arguments:
%       inDir : starting directory, e.g. '' or 'mydir' or 'mydir/'
%       fNames: nx1 cell array to receive names
%       pars.displayOn set to 1 to print out progress
%       pars.extensions cell array of extensions, e.g. {'.tif' '.mrc'}
%       pars.maxDepth: how far below initial directory level.

searchDir=AddSlash(inDir);
if numel(inDir)==0
    searchDir='./';
end;
        if pars.displayOn
            disp(['Searching ' searchDir]);
        end;
d=dir(searchDir);
nf0=size(fNames,1);
nf=nf0;
for i=1:numel(d)
    if d(i).name(1)=='.' % ignore ridiculous names
        continue;
    end;%                 disp(d(i).name)
    if d(i).isdir && pars.maxDepth>0    % new directory; step into it
        inDirNext=AddSlash([searchDir d(i).name]);
        pars1=pars;
        pars1.maxDepth=pars.maxDepth-1;
        fNames=SearchDirectories(inDirNext,fNames,pars1);
    else  % it's a file
        name=d(i).name;
%         disp(name);
        for j=1:numel(pars.extensions)  % check for a valid extension
%            numC=numel(pars.extensions{j});
%             ok=numel(name)>=numC && strndcmp(name,pars.extensions{j},numC);
                        ok=strndcmp(name,pars.extensions{j});
            if ok
                break;
            end;
        end;
        if ok
            nf=nf+1;
            fNames{nf,1}=[inDir name];
        end;
    end; % if isdir
end % for
if pars.displayOn && nf>nf0
    disp(['  ' num2str(nf-nf0) ' file' ess(nf-nf0) ' found in ' inDir]);
end;
return
