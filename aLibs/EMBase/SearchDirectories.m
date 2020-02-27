function fNames=SearchDirectories(inDir,fNames,pars)
% **** There is a newer version on farnam 10/11/19
% Recursive search for files having the correct endings or
% regexp pattern. Builds up the
% cell array of complete filenames fNames(:,1) and the directory entry structs fNames(:,2).
% Arguments:
%       inDir : starting directory, e.g. '' or 'mydir' or 'mydir/'
%       fNames: cell array (can be empty) to which names are added
%       pars.displayOn set to 1 to print out progress
%       pars.extensions cell array of extensions, e.g. {'.tif' '.mrc'}
%         or, alternatively,
%       pars.regexps a cell array of reg expression strings, in descending
%         order.
%       pars.maxDepth: how far below initial directory level to search

if nargin<3
    pars=struct;
end;
if nargin<2
    fNames={};
end;
defPars.displayOn=0;
defPars.extensions={'.m'};
defPars.regexps={};
defPars.maxDepth=4;
pars=SetDefaultValues(defPars,pars);

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
    elseif ~d(i).isdir  % it's a file
        name=d(i).name;
        %         disp(name);
%         if isfield(pars.regexps{pars.maxDepth}) && numel(pars.regexps)>0 
        if isfield(pars,'regexps') && pars.maxDepth>0 && numel(pars.regexps)>=pars.maxDepth
                % use regular expression
            ok=numel(regexp(name,pars.regexps{pars.maxDepth}))>0;
        else
            for j=1:numel(pars.extensions)  % check for a valid extension
                %            numC=numel(pars.extensions{j});
                %             ok=numel(name)>=numC && strndcmp(name,pars.extensions{j},numC);
                ok=strndcmp(name,pars.extensions{j});
                if ok
                    break;
                end;
            end;
        end;
        
        if ok
            nf=nf+1;
            fNames{nf,1}=[inDir name];
            fNames{nf,2}=d(i); % pick up the whole struct
        end;
    end; % if isdir
end % for
if pars.displayOn && nf>nf0
    disp(['  ' num2str(nf-nf0) ' file' ess(nf-nf0) ' found in ' inDir]);
end;
return
