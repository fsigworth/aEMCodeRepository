function ExpandZMovies(inDir,outDir)
% use:
%  ExpandZMovies   (put up file selectors)
%  ExpandZMovies('myInputDirectory/','myOutputDirectory/');
% 
% Recursively search a directory tree, picking all files having the proper
% extension.  We expand them and write them into a parallel tree (or the
% same tree).  An input file namez.tif gets written as name.mrcs (or
% whatever is given for outExtension).  Alternatively, if useOriginalName
% is set, the original name is read from the ZTiff file and used.
% If no arguments are given the program puts up file selectors; however
% in batch mode if inDir and outDir aren't given, they are simply set to 
% the current directory.

useParFor=0;  % set to run parallel processes
overwrite=1;  % replace existing output files
inExtension='z.tif';  % search for filename ending with this
outExtension='.mrcs';  % replace with this extension.
useOriginalName=0;  % If true, use the original name stored in the file.
batchMode=0;

% Force batch mode if we have no GUI
batchMode=batchMode || ~usejava('desktop');

classes={'uint8' 'int16' 'single' '' '' '' 'uint16'};
if nargin<1
    if batchMode
        inDir='./';
    else
        inDir=uigetdir('','Select the input directory');
        if isnumeric(inDir)
            return
        end;
    end;
    if nargin<2 && batchMode
        outDir=inDir;
    else
        outDir=uigetdir(inDir,'Select the output directory');
        if isnumeric(outDir)
            return
        end;
    end;
end;
inDir=AddSlash(inDir);
outDir=AddSlash(outDir);

fileNames=SearchDirectories(inDir,outDir,{});

nfiles=size(fileNames,1);
disp(['There are ' num2str(nfiles) ' files to expand']);

if useParFor
    tic
    parfor fi=1:nfiles
        finName=fileNames{fi,1};
        foutName=fileNames{fi,2};
        foutPath=fileNames{fi,3};
        [m,s,fok]=ReadZTiff(finName);
        if fok
            if ~exist(foutPath,'dir')
                disp(['Creating the directory ' foutPath]);
                mkdir(foutPath);
            end;
            if useOriginalName && isfield(s,'origName') && numel(s.origName)>1
                foutName=[foutPath s.origName];
            end;
            
            if ~isfield(s,'origClass') || numel(s.origClass)<1
                s.origClass='single';
            end;
            classOk=strcmp(classes,s.origClass);
            if any(classOk)
                outMode=find(classOk,1)-1;
            else
                outMode=2;  % default is single
            end;
            
            fprintf(1,'Writing %s\n',foutName);
            WriteMRC(m,s.pixA,foutName,outMode);
        end;
    end;
    toc
else
    tic
    for fi=1:nfiles
        finName=fileNames{fi,1};
        foutName=fileNames{fi,2};
        foutPath=fileNames{fi,3};
        [m,s,fok]=ReadZTiff(finName);
        if fok
            if ~exist(foutPath,'dir')
                disp(['Creating the directory ' foutPath]);
                mkdir(foutPath);
            end;
            if useOriginalName && isfield(s,'origName') && numel(s.origName)>1
                foutName=[foutPath s.origName];
            end;
            
            if ~isfield(s,'origClass') || numel(s.origClass)<1
                s.origClass='single';
            end;
            classOk=strcmp(classes,s.origClass);
            if any(classOk)
                outMode=find(classOk,1)-1;
            else
                outMode=2;  % default is single
            end;
            
            fprintf(1,'Writing %s\n',foutName);
            WriteMRC(m,s.pixA,foutName,outMode);
        end;
    end;
    toc
end;

if isdeployed && pars.displayOn
    pause(2);
    close;  % close the figure window.
end;

    function fNames=SearchDirectories(inDir,outDir,fNames)
        % Recursive search for files having the correct extensions.  Builds up the
        % list of complete filename fNames.
        d=dir(inDir);
        for i=1:numel(d)
            if d(i).name(1)~='.' % ignore ridiculous names
                if d(i).isdir    % new directory; step into it
                    inDirNext=AddSlash([inDir d(i).name]);
                    outDirNext=AddSlash([outDir d(i).name]);
                    disp(['Searching ' inDirNext]);
                    fNames=SearchDirectories(inDirNext,outDirNext,fNames);
                else  % it's a file
                    name=d(i).name;
                    if strndcmp(name,inExtension);
                        outName=[name(1:end-numel(inExtension)) outExtension];
                        if overwrite || ~exist([outDir outName],'file')
                            nf=size(fNames,1)+1;
                            fNames{nf,1}=[inDir name];
                            fNames{nf,2}=[outDir outName];
                            fNames{nf,3}=outDir;
                        else
                            disp(['  already expanded: ' outDir outName])
                        end;
                    end;
                    
                end;
            end;
        end;
        
    end
end
