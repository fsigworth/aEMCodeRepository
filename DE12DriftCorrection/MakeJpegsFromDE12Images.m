% MakeJpegsFromDE12Images.m
% Accumulate DE camera raw images, and perform dark and gain correction.
% We write the files AccumSum.mat and AccumSumOdd into the raw image folder,
% along with jpeg and mrc corrected images in the enclosing folder.
%
% For computing our own dark and gain references we assume
% a folder structure such as
%   expt/Bright
%   expt/Dark
%   expt/Images (can be any name, contains actual data images)
% 
% Each of these contains the DE camera raw-image folders, so for example
% the contents of the Bright folder might be
%   expt/Bright/DE_20120113_001135_358/         -- one movie directory
%   expt/Bright/DE_20120113_001135_358/RawImage_0.tif
%   expt/Bright/DE_20120113_001135_358/RawImage_1.tif
% and more raw images.  There are typically more movie directories as well.
% 
% At the end we write out the following files:
%   expt/Bright/DE_20120113_001135_358/isum.mat -raw sum of images in movie
%   expt/Bright/AccumSum.mat      --raw sum of all movies in expt/Bright/
%   expt/Bright/AccumSumOdd.mat   --raw sum of every other movie
% 
% Once we have done this for the Bright and Dark folders, their AccumSums
% are used to do correction of each movie in the Images folder.  An example
% of files written there for one movie:
%   expt/Images/DE_20111229_231551_002/isum.mat
%   expt/Images/DE_20111229_231551_002.mrc
%   expt/Images/DE_20111229_231551_002.jpg
% 
% Matlab will follow symlinks but not MacOS aliases, for example to pick up
% Bright and Dark folders; e.g.
%    ln -s favorites/Bright expt/Bright
% 
% Note that MRC files are written after operating on the image m with
% rot90(m,3).  This puts m into Cartesian coordinates m(x,y) to be
% consistent with the EM standards.

useExistingAccums=1;  % Don't create new files if old ones exist.
writeMRC=1;    % Write out an MRC of the corrected image too.
epsi=1e-3;  % Fraction of histogram that is truncated for jpeg display
binFactor=2; % Binning for jpeg file
skipFrames=1; % Number of frames to skip at start of movie
pixA=2.9;  %Written into the MRC files.

% Where to find the references if we don't have 'bright' and 'dark' folders
drefPath='~/EMWork/Hideki/111229/111229/AMPAR_boxB_slot3/ReferenceBuffers';
darkRefName='DE_20111227_005853_161_DarkFrame.tif';
grefPath='~/EMWork/Hideki/111229/111229/AMPAR_boxB_slot3/ReferenceBuffers';
gainRefName='DE_20111227_010026_059_GainFrame.tif';


% First, we have the user select the enclosing folder (expt/ in the example
% above.
uiPath=pwd;
uiPath=uigetdir(uiPath,'Select a folder containing Images/');
if isa(uiPath,'numeric')
    return
end;
cd(uiPath);
folders=FindFilenames('.','\w+',1);  % Get all directories with any names
origFolders=folders;  % copy for debugging.
char(folders)

% txt=inputdlg('Pixel size in Å');  % Get the pixel size from a dialog box.
% if numel(txt)>0
%     pixA=str2num(char(txt));
% else
%     pixA=0;
% end;

% pixA=input('Pixel size in Å? ');  % Get it from the command window.

figure(1);
SetGrayscale;
subplot(2,1,1);

% See if we have a dark reference folder
drn=0;
q=strcmpi(folders,'dark');
if any(q)  % There is a 'dark' folder;
    ind=find(q,1,'first');
    dname=folders{ind};  % get its name
    folders(ind)=[];  % remove that name from the folders list.
    % Look into the 'dark' folder for 'AccumSum.mat'
    r=FindFilenames(dname,'(?i)AccumSum\.mat');
    if numel(r)>0 && useExistingAccums
        accumName=[AddSlash(dname) r{1}];
        disp(['Loading dark sum: ' accumName]);
        load(accumName);  % get mAccum, nAccum
    else  % no AccumSum file: try to make one.
        disp('Accumulating dark images');
        [mAccum, nAccum]=deGetAccumSum(dname,skipFrames);
    end;
    drm=mAccum;
    drn=nAccum;
end;
if drn>0
    dref=drm/drn;
else  % read the default
    disp(' using the default dark reference');
    dref=single(imread([AddSlash(drefPath) darkRefName]));
end;

% See if we have a bright reference folder
brn=0;
q=strcmpi(folders,'bright');
if any(q)  % There is a 'bright' folder;
    ind=find(q,1,'first');  % get its name
    bname=folders{ind};
    folders(ind)=[];  % remove that name from the folders list.
    % Look into the 'bright' folder for 'AccumSum.mat'
    r=FindFilenames(bname,'(?i)AccumSum\.mat');
    if numel(r)>0 && useExistingAccums
        accumName=[AddSlash(bname) r{1}];
        disp(['Loading bright sum: ' accumName]);
        load(accumName);  % get mAccum, nAccum
    else  % no AccumSum file: try to make one.
        disp('Accumulating bright images');
        [mAccum, nAccum]=deGetAccumSum(bname,skipFrames);
    end;
    brm=mAccum;
    brn=nAccum;
end;
if brn>0
    gref=brm/brn-dref;
    gref=gref/mean(gref(:));
else  % read the default
    disp(' using the default gain reference');
    gref=single(imread([AddSlash(grefPath) gainRefName]));
end;

% Go through all the remaining folders
if numel(folders)>0
    disp('Operating on data:');
end;
for j=1:numel(folders)
    md=AddSlash(folders{j});  % A folder containing movie folders
    r=FindFilenames(md,'DE',1);  % Get movie folder names starting with DE
    for i=1:numel(r)  % Operate on one movie folder
        movieDir=[md AddSlash(r{i})];
        disp(movieDir);
        subplot(2,1,1);
        [m nim]=deGetRawSum(movieDir,skipFrames);
        if nim>0
            corrImage=(m-nim*dref)./gref;
            outName=[md r{i}];
            subplot(2,1,2);
            imagesc(corrImage);
            title([outName '.mrc'],'interpreter','none');
            drawnow;
            if writeMRC
%                 WriteMRC(rot90(corrImage,3),pixA,[outName '.mrc']);
                WriteMRC((corrImage)',pixA,[outName '.mrc']);
            end;
            WriteJpeg(corrImage,[outName '.jpg']);
        end;
    end;
end;

