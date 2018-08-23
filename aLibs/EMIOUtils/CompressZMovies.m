function CompressZMovies(inDir,outDir,pars)
% Recursively search a directory tree, picking all files having the proper
% extensions.  We compress them and write them into a parallel tree (or the
% same tree).  An input file name.mrcs gets written as namez.tif
% If no arguments are given, the program puts up file selectors.

lp=struct;
lp.useParFor=0;
lp.overwrite=1;
lp.extensions={'.mrc' '.mrcs'};
lp.padK2Movies=0;
lp.snrRatio=200;  % target noise-to-noise ratio
lp.lfCutoff=.05; % radius of frequency domain not fitted
lp.displayOn=0;
lp.falconMasking=0;
% For merged images:
%pars.snrRatio=200;  % target noise-to-noise ratio
%pars.lfCutoff=.2; % radius of frequency domain not fitted
if nargin<3
    pars=struct;
end;
pars=SetOptionValues(lp,pars);


if nargin<2
    inDir=uigetdir('','Select the input directory');
    if isnumeric(inDir)
        return
     end;
    outDir=uigetdir(inDir,'Select the output directory');
    if isnumeric(outDir)
        return
    end;
    pars.displayOn=(isdeployed || pars.displayOn) && ~useParFor;  % we put up the figure window if running alone.
end;
inDir=AddSlash(inDir);
outDir=AddSlash(outDir);

fileNames=SearchDirectories(inDir,outDir,cell(0,1));

nfiles=size(fileNames,1);
disp(['There are ' num2str(nfiles) ' files to compress']);

if pars.useParFor
    tic
    parfor fi=1:nfiles
        lpars=pars;
        lpars.displayOn=0;
        finName=fileNames{fi,1};
        foutName=fileNames{fi,2};
        [m,s,fok]=ReadMovie(finName);
%         if doFixK2
            
        if fok
            pa=fileparts(foutName);
            if ~exist(pa,'dir')
                disp(['Creating the directory ' pa]);
                mkdir(pa);
            end;                
            fprintf(1,'Writing %s\n',foutName);
            [~,nm,ex]=fileparts(finName);  % we'll strip the path
            lpars.origName=[nm ex];
            WriteZTiff(m,s.pixA,foutName,lpars);
        end;
    end;
    toc
else
    for fi=1:nfiles
        finName=fileNames{fi,1};
        foutName=fileNames{fi,2};
        tic
        [m,s,fok]=ReadMovie(finName);
        if pars.falconMasking
            m=f2FixEdges(m);
        end;
        if pars.padK2Movies
            me=mean(m(:));
            n1=NextNiceNumber(size(m,1));
            m=Crop(m,n1,1,me);
        end;
        if fok
            pa=fileparts(foutName);
            if ~exist(pa,'dir')
                disp(['Creating the directory ' pa]);
                mkdir(pa);
            end;                
            fprintf(1,['Writing ' foutName '...']);
            [~,nm,ex]=fileparts(finName);  % we'll strip the path
            pars.origName=[nm ex];
            WriteZTiff(m,s.pixA,foutName,pars);
            fprintf(1,' %gs\n',toc);
        end;
    end;
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
                    for j=1:numel(pars.extensions)  % check for a valid extension
                        numC=numel(pars.extensions{j});
                        ok=numel(name)>numC && strcmpi(name(end-numC+1:end),pars.extensions{j});
                        if ok
                            break;
                        end;
                    end;
                    if ok
                        outName=[name(1:end-numC) 'z.tif'];
                        if pars.overwrite || ~exist([outDir outName],'file')
                            nf=size(fNames,1)+1;
                            fNames{nf,1}=[inDir name];
                            fNames{nf,2}=[outDir outName];
                        else
                            disp(['  already compressed: ' outDir outName])
                        end;
                    end;
                end;
            end;
        end;
        
    end
end
