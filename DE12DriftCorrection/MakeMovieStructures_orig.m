% deMakeMovieStructures.m
% Load each DE camera movie found in expt/Images, where expt is the base
% directory for the experiment.  Create for it the s
% data structure which is stored in a folder expt/Tracking
% s has the fields
%
%   s.groupN % number of raw images summed to make a group or 'frame'
%   s.dr1  % two independent dark images, scale is the same as a frame
%   s.dr2
%   s.br1  % two indpendent averaged bright images, scaled like a frame
%   s.br2
%   s.imgs  % stack of frames, uint16 values, uncorrected
%   s.frameDose  % e/A^2 in each frame
%
% To make a fully corrected image if frame i using reference group1,
%   m=(single(s.imgs(:,:,i))-s.dr1)./(s.br1-s.dr1)*mean(s.br1(:));

doWriteFiles=1; % Write out the s structures
% pixA is set here
pixA=2.9;

groupN=5;  % must be at least 2 if we want ramp compensation
% There must be at least 3 groups for processing, since we sum groups
% starting with number 2 to compute means.

refScale=groupN;
maxNgrps=60;  % don't process more than this.
divisor=1;
skipN=1;     % number of raw frames to skip
stride=1;   % successive raw frame stride
istep=1;    % increment through directory

rexp='DE_20.+';  % raw frame directory name is DE_20 followed by one or more characters.
deadColumns=[];
refBaseName='';
dirStart=1;
outPostscr='';

figure(1);
clf;
SetGrayscale;
doFlipSum=0;
doFlipRefs=0;

mainPath=AddSlash(uigetdir());
if numel(mainPath)<3
    return
end;
cd(mainPath);
mrcRefs=0;
dpath='Dark/';
bpath='Bright/';
refPath=mainPath;
imagePath='Images/';

disp(['Image path: ' imagePath]);

% Where to put -s files.  Create a directory if needed.
outPath=[mainPath 'Tracking/'];
if doWriteFiles && ~DirectoryExists(outPath)
    mkdir(outPath);
    disp(['Directory was created: ' outPath]);
end;

dirNames=FindFilenames(imagePath,rexp,1);  % search for directories
disp([num2str(numel(dirNames)) ' movies found']);

% Get references and create the base structure.  The raw accumulated images
% are rotated by 270 degrees, as are the raw images themselves.

s0=struct;
s0.groupN=groupN;
s0.pixA=pixA;

dref=load([refPath dpath 'AccumSum.mat']);
drefo=load([refPath dpath 'AccumSumOdd.mat']);

s0.dr1=rot90((dref.mAccum-drefo.mAccumOdd)/(dref.nAccum-drefo.nAccumOdd)*groupN,3);
s0.dr2=rot90(drefo.mAccumOdd/drefo.nAccumOdd*groupN,3);

bref=load([refPath bpath 'AccumSum.mat']);
brefo=load([refPath bpath 'AccumSumOdd.mat']);
s0.br1=rot90((bref.mAccum-brefo.mAccumOdd)/(bref.nAccum-brefo.nAccumOdd)*groupN,3);
s0.br2=rot90(brefo.mAccumOdd/brefo.nAccumOdd*groupN,3);

s0.dr1=deFixColumns(s0.dr1,deadColumns);
s0.dr2=deFixColumns(s0.dr2,deadColumns);
s0.br1=deFixColumns(s0.br1,deadColumns);
s0.br2=deFixColumns(s0.br2,deadColumns);

disp('References loaded');

infoSearchString='(electron/pixel/frame)=';

for idir=dirStart:istep:numel(dirNames)
    imName=dirNames{idir};
    imPath=AddSlash(imName);
    
    s=s0;
    
    % Note that in each read operation we rotate by 270 degrees.
    
    %     [ny nx]=size(dref.mAccum);
    
    %%
    disp(['reading ' imPath]);
    
    % pick up the dose from the info.txt file.
    pixelDose=0;
    infoHandle=fopen([imagePath imPath 'info.txt']);
    if infoHandle>0  % file exists
        lines=textscan(infoHandle,'%s','delimiter','\n');
        fclose(infoHandle);
        for j=1:numel(lines{1})
            line=lines{1}{j};
            q=strfind(line,infoSearchString);
            if numel(q)>0
                pixelDose=str2num(line(q(1)+numel(infoSearchString):numel(line)));
                break
            end;
        end;
    end;
    frameDose=pixelDose*groupN/(pixA^2)
    s.frameDose=frameDose;
    
    names=deGetRawImageNames([imagePath imPath]);
    nim=numel(names)-skipN;
    ngrps=floor(nim/(groupN*stride));
    ngrps=min(maxNgrps,ngrps);
    % Read the first frame to get the image size
    m1=imread(names{1});
    [ny nx]=size(m1);
    
    s.imgs=uint16(zeros(nx,ny,ngrps));
    
    for i=1:ngrps;
        m=uint16(zeros(nx,ny));
        for j=1:stride:groupN*stride
            ind=(i-1)*groupN*stride+j+skipN;
            name=names{ind};
            m1=uint16(ReadEMFile(name));
            m=m+m1;
        end;
        if doFlipSum
            m=flipud(m);
        end;
        %         m=m-mean(m(:));
        % We fix the dead columns and then convert back to integer.
        s.imgs(:,:,i)=uint16(deFixColumns(single(m),deadColumns));
        mcorr=(single(m)-s.dr1)./(s.br1-s.dr1);
        imac(uint8(imscale(BinImage(mcorr,4),256,1e-4)));
        title(i);
        drawnow;
    end;
    %%
    
    outName=[outPath imName outPostscr];
    if doWriteFiles
        %%
        disp(['Writing ' outName]);
        drawnow;
        save([outName '-s.mat'],'s');
    end;
end;


