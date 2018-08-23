function fNames=SearchDirectories(inDir,fNames,pars)
        % Recursive search for files having the correct extensions.  Builds up the
        % cell array of complete filenames fNames.
%         arguments:
%       inDir : starting directory, e.g. '.' but not ''
%       fNames: nx1 cell array to receive names
%       pars.displayOn set to 1 to print out progress
%       pars.extensions cell array of extensions, e.g. {'.tif' '.mrc'}
        d=dir(inDir);
        for i=1:numel(d)
            if d(i).name(1)~='.' % ignore ridiculous names
                if d(i).isdir    % new directory; step into it
                    inDirNext=AddSlash([inDir d(i).name]);
                    if pars.displayOn
                        disp(['Searching ' inDirNext]);
                    end;
                    fNames=SearchDirectories(inDirNext,fNames,pars);
                else  % it's a file
                    name=d(i).name;
                    for j=1:numel(pars.extensions)  % check for a valid extension
                        numC=numel(pars.extensions{j});
                        ok=numel(name)>numC && strcmpi(name(end-numC+1:end),pars.extensions{j});
                        if ok
                            break;
                        end;
                    end;
                    if ok
                            nf=size(fNames,1)+1;
                            fNames{nf,1}=[inDir name];
                    end;
                end; % if isdir
            end;  % if ~='.'
        end % for
end
