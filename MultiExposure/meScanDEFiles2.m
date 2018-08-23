function infos=meScanDEFiles2(mi0,numExposures)
% Look at the movie-file directories in the given directory (concatenation of basePath and imagePath)
% and create a micrograph info structure mi for each set of exposures.
% mi.imagePath
% mi.procPath
% mi.infoPath
% mi.baseFilename
% Save the structure in a file having the full name
% [infodir basename 'minfo.mat']
% If the file already exists, do not overwrite it but instead load its
% contents as mi.
% Then return a cell array of all the structures.

d=dir([mi0.basePath mi0.imagePath]);

nds=numel(d);
%  A typical filename is 'DE_20111229_233334_160.mrc'.
%  The regexp pattern is
pattern='DE_\d+_\d+_.+\.mrc';

% Scan the directory and pick up the names
infos={};
ind=0;
j=0;
ok=1;
while (ind < nds) && ok
    %     Try to pick up a set of filenames
    [ind ok]=GetNextImageFilename(d,ind,pattern);  % Get first image
    if ~ok
        break  % immediately exit the loop
    end;
    name=d(ind).name;
    mi=mi0;
    [pa baseFilename]=fileparts(name);
    mi.baseFilename=baseFilename;
    mi.imageFilenames{1}=name;
    for i=2:numExposures
        [ind ok]=GetNextImageFilename(d,ind,pattern);
        if ~ok
            break  % exit for loop
        end;
        mi.imageFilenames{i}=d(ind).name;
    end;
    if ok  % we successfully got the right number of files
        j=j+1;
        infos{j}=mi;
    end;
end % while
end

    function [ind ok]=GetNextImageFilename(d,ind,pattern)
        ind=ind+1;  % move ahead
        ok=1;
        nd=numel(d);
        while ind<=nd
            name=d(ind).name;
            if ~d(ind).isdir && numel(regexp(name,pattern))
                return
            else
                ind=ind+1;
            end;
        end;
        ok=ind<nd;
    end

