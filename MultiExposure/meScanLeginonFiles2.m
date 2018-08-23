function infos=meScanLeginonFiles2(imagedir,procdir, infodir, overwrite)
% function infos=meScanLeginonFiles(serverdir,localdir,infodir)
% Look at the files in the given directory (concatenation of serverdir and
% localdir) and create a micrograph info structure mi for each set of exposure
% files.  The following fields are filled in:
% mi.imagePath
% mi.procPath
% mi.infoPath
% mi.baseFilename
% Save the structure in a file having the full name
% [infodir basename 'minfo.mat']
% If the file already exists, do not overwrite it but instead load its
% contents as mi.
% Then return a cell array of all the structures.

if nargin<4
    overwrite=0;
end;

% Create the maps for exposure types
exTypes=  {'en','em','ef'};
exSequence=containers.Map(exTypes,[ 1 2 3 ]);



% serverdir='/Volumes/TetraData/EMWork/Liguo/FavoriteBK/';  % main directory where images are stored
% % serverDir='/Volumes/fred/EMWork/Liguo/';  % main directory where images are stored
% imagedir='10dec18a/';         % subdirectory for these images
% infodir='10dec18aMi/';            % complete path to where the info files are stored.
% procdir=imagedir;
% overwrite=1;

% scan the current directory, looking for Leginon exposure files.  For each
% set of exposures, create a micrograph info structure, store it as an
% XXXminf.mat file, and also add it to the infos array of structures.

d=dir(imagedir);
nds=numel(d);
% We identify valid image file names as ones having the pattern
%   *[0-9]e[A-Za-z]*.mrc
% So that an example image name could end with '0001ef.mrc'
% We then identify all related files as those having the same base, i.e.
% having the same name up to the final numeric character.
UniqueBases={};
Bases=cell(nds,1);
% ImTypes=cell(nds,1);
Names=cell(nds,1);
Exts=cell(nds,1);
ImSeq=zeros(nds,1);
% Scan the directory and pick up the names and exposure type (seq) and
% extensions
for i=3:nds
    [pa base imtype ext]=LeginonFileParts(d(i).name);
    %     base
    %     imtype
    Names{i}=d(i).name;
    Bases{i}=base;
    if isKey(exSequence,imtype)
        seq=exSequence(imtype);
    else
        seq=0;
    end;
    ImSeq(i)=seq;
    Exts{i}=ext;
    if strcmp(ext,'.mrc') && seq>0  % an exposure file
        if ~any(strcmp(base, UniqueBases))  % not the same as the previous base name
            UniqueBases=[UniqueBases ; base];
        end;
    end;
end;
nb=numel(UniqueBases);

mi=meCreateMicrographInfoStruct10;

% Now scan for defocus pairs or triples, sort their names, and insert them
% into the mi structure.
infos={};
for i=1:nb
    infoFilename=[infodir UniqueBases{i} 'mi.mat'];
%     fi=fopen(infoFilename,'r');
    %     if fi>0 && ~overwrite % the file exists already
    if FileExists(infoFilename) && ~overwrite % use the existing file
        disp([infoFilename ' -exists']);
        load(infoFilename);  % load the existing one.
    else    % construct the mi struct from scratch
        % First, find all the files with a given base name
        q=strcmp(UniqueBases{i},Bases);
        if any(q)
            %     relfiles={Names{q}};
            relatedFiles={Names{q}};
            %             relatedFileTypes={ImTypes{q}};
            relatedFileSeq=ImSeq(q);
            relatedFileExts={Exts{q}};
            
            % Fill fields in the micrograph info structure
            %             mi=meCreateMicrographInfoStruct;
            mi.baseFilename=UniqueBases{i};
            mi.imagePath=imagedir;
            mi.procPath=procdir;
            mi.infoPath=infodir;
            
            % Find out how many exposures we have, including extension '.mrc'
            mrcs=strcmp('.mrc',relatedFileExts);
            q=find(mrcs & relatedFileSeq');
            relatedImages={relatedFiles{q}};
            imageSeq=relatedFileSeq(q);
            [vals seq]=sort(imageSeq);
            nExposures=numel(seq);
            if nExposures<1
                error(['code error: no valid exposure files found for '...
                    mi.baseFilename]);
            end;
            mi.imageFilenames={relatedImages{seq}};
            %             char(mi.rawFile)
            %
            % %
            %             switch nExposures
            %                 case 0
            %                     error('code error: no valid exposures files found.')
            %                 case 1
            %                     mi.rawFile=relatedFiles{exposureIndices};
            %                 case 2
            %                     mi.rawFile={{[mi.baseFilename 'en.mrc']};...
            %                         {[mi.baseFilename 'ef.mrc']}};
            %                 case 3
            %                     mi.rawFile={{[mi.baseFilename 'en.mrc']};...
            %                         {[mi.baseFilename 'em.mrc']};...
            %                         {[mi.baseFilename 'ef.mrc']}};
            %                 otherwise
            %                     error(['Too many exposure files found in ' mi.baseFileName]);
            %             end;
            save(infoFilename,'mi');
            disp([infoFilename ' -created'])
        end;
    end;
    infos{i}=mi;
end;

