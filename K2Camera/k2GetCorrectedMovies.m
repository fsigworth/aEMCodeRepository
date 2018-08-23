% k2GetCorrectedMovies.m
% Read movie files as in k2DriftTracker but write out gain-corrected movies
% padded and rotated the same as the summed images produced by
% k2DriftTracker.
nrot=3;  % number of 90 degree rotations to match the gain reference image

[names, pathName]=uigetfile('*mi.*','Select mi files','multiselect','on');
if isnumeric(pathName) % File selection was cancelled
    return
end;
if ~iscell(names)  % if only one is selected, it's a string not a cell array
    names={names};
end;
[rootPath, infoPath]=ParsePath(pathName);
cd(rootPath);

for miIndex=1:numel(names)
    miName=[infoPath names{miIndex}];
    disp(['Reading ' miName]);
    mi=ReadMiFile(miName);
    
    flagSegments=mi.frameSets;
    numDefs=max(1,size(flagSegments,1));
    
    gainRef=ReadEMFile(mi.gainRefName);
    movieName=[mi.moviePath mi.movieFilename];
    [pa,nm]=fileparts(mi.movieFilename);
    for segIndex=1:numDefs  % defocus stretch
        disp(['Reading ' movieName ' (segment ' num2str(segIndex) ')']);
        %     Get the header information
        [mv, s]=ReadMovie(movieName,1,0);
        n=[1 1]*NextNiceNumber(max(s.nx,s.ny)); %  final frame size =[3840 3840] for k2 data
        inds=flagSegments(segIndex,:); % of the form [startIndex endIndex]
        inds=min(s.nz,inds);
        mi.frameSets(segIndex,:)=inds;  % trim the frame numbers.
        nim=inds(2)-inds(1)+1;
        disp([' frames ' num2str(inds)]);
        if nim < 1
            error('no frames in this segment');
        else
            m=single(rot90(ReadMovie(movieName,inds(1),nim),nrot));
            m0=m.*repmat(gainRef,1,1,nim);
        end;
        mi.imageSize=n;
        m1=Crop(m0,n,1,mean(m0(:))); % pad the image with mean-value pixels
        corrMovieName=[mi.moviePath nm '_corr_seg_' num2str(segIndex) '.mrc'];
        disp(['Writing ' corrMovieName]);
        WriteMRC(m1,mi.pixA,corrMovieName);
    end;
end;
