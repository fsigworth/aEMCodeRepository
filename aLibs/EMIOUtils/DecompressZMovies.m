function DecompressZMovies
% Recursively search a directory tree, picking all files having the z.tif
% suffix.  We decompress them and write them into a parallel tree (or the
% same tree).  An input file namez.tif gets written as name.mrc.  The mode
% of the mrc file is 2 (=single) by default, but can be set to the origMode
% of the compressed data.

useParFor=0;
overwrite=1;
extension='z.tif';
outExtension='.mrc';
forceMode=false;  % force the output to be mode 2 (single precision)

inDir=uigetdir('','Select the input directory');
if isnumeric(inDir)
    return
else
    inDir=AddSlash(inDir);
end;
outDir=uigetdir(inDir,'Select the output directory');
if isnumeric(outDir)
    return
else
    outDir=AddSlash(outDir);
end;

fileNames=SearchDirectories(inDir,outDir,{});

nfiles=size(fileNames,1);
disp(['There are ' num2str(nfiles) ' files to decompress']);

if useParFor
    tic
    parfor fi=1:nfiles
        finName=fileNames{fi,1};
        foutName=fileNames{fi,2};
        [m,s,fok]=ReadZTiff(finName);
        mode=2;
        if isfield(s,'origClass') && ~forceMode
            switch s.origClass
                case 'uint8'
                    mode=0;
                case 'int16'
                    mode=1;
            end;
        end;
        if fok
            fprintf(1,'Writing %s\n',foutName);
            WriteMRC(m,s.pixA,foutName,mode);
        end;
    end;
    toc
else
    for fi=1:nfiles
        finName=fileNames{fi,1};
        foutName=fileNames{fi,2};
        tic
        [m,s,fok]=ReadZTiff(finName);
        mode=2;
        if isfield(s,'origClass') && ~forceMode
            switch s.origClass
                case 'uint8'
                    mode=0;
                case 'int16'
                    mode=1;
            end;
        else
            mode=2;
        end;
        if fok
            fprintf(1,['Writing ' foutName '...']);
            WriteMRC(m,s.pixA,foutName,mode);
            fprintf(1,' %gs\n',toc);
        end;
    end;
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
                    if ~exist(outDirNext,'dir')
                        disp(['Creating ' outDirNext]);
                        mkdir(outDirNext);
                    else
                        disp(['Searching ' outDirNext]);
                    end;
                    fNames=SearchDirectories(inDirNext,outDirNext,fNames);
                else  % it's a file
                    name=d(i).name;
                    el=numel(extension);
                    if numel(name)>el && strcmpi(name(end-el+1:end),extension);
                        outName=[name(1:end-el) outExtension];
                        if overwrite || ~exist([outDir outName],'file')
                            nf=size(fNames,1)+1;
                            fNames{nf,1}=[inDir name];
                            fNames{nf,2}=[outDir outName];
                        end;
                    end;
                end;
            end;
        end;
    end
end

