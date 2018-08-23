function [si,imgs,fnames,siPath,altImgs]=reLoadStackFiles(opts,fnames,siPath)
% Given filenames, read the *si.mat files and merge the si structures
% together.  Also load corresponding stacks *stack.mrc' and concatenate
% them into imgs.  If desired, the corresponding unsubtracted stacks are
% also concatenated to form the altImgs.
% Both input arguments are optional. fnames can be either a
% string or a cell array of strings, and if absent or empty, a file selector is put
% up.  If no file was selected, si is returned as an empty struct and imgs
% is an empty matrix.

if nargin<1
    opts=struct;
end;
if nargin<2
    fnames=[];
end;
if nargin<3
    siPath='';
end;
readAltImgs=(nargout>4);

if ~isfield(opts,'doNormalize')
    opts.doNormalize=0;
end;

si=struct;
imgs=[];

if numel(fnames)<1
    % Put up a file selector for files *si.mat,
    disp('Getting si file names');
    [fnames, siPath]=uigetfile('*si.mat','Select si files','multiselect','on');
    if isnumeric(siPath)  % user clicked Cancel
        return
    end;
    pa=ParsePath(siPath);  % go to main directory
    cd(pa);
    siPath=AddSlash(siPath);
end;
if ~iscell(fnames)
    fnames={fnames};
end;

%%
% read and merge the
% structures into si.  Also load corresponding stacks *stack.mrc' and concatenate
% them into imgs.
allImgs=[];
allSi=struct;
allAltImgs=[];
for i=1:numel(fnames)
    siName=fnames{i};
    disp(['Loading ' siName]);
    p=regexp(siName,'si.mat');
    if numel(p)<1
        error(['No match found in name ' siName]);
    end;
%     si=LoadStruct([siPath siName]);
    load([siPath siName]);
    stackName=[siName(1:p(end)-1) 'stack.mrc'];
    disp(stackName);
    imgs=ReadEMFile([siPath stackName]);
    [allSi, allImgs]=rsStackConcatenate(si,imgs,allSi,allImgs);
    if readAltImgs
        altStackName=[siName(1:p(end)-1) 'ustack.mrc'];
        disp(altStackName);
        altImgs=ReadEMFile([siPath altStackName]);
        %         Concatenate images only
        szAll=min(size(allAltImgs,3),numel(allAltImgs));
        sz1=size(altImgs,3);
        if szAll<1
            allAltImgs=altImgs;
        else
            allAltImgs(:,:,szAll+1:szAll+sz1)=altImgs;
        end;
    end;
end;

if ~isfield(allSi,'mbnOffset')
    allSi.mbnOffset=defaultMbnOffsetA/allSi.pixA;
    disp(['Using the default mbn offset: ' num2str(allSi.mbnOffset*allSi.pixA) ' A']);
end;
si=allSi;

n0=size(allImgs,1);
nImgs=size(allImgs,3);

imgs=allImgs;
altImgs=allAltImgs;

if opts.doNormalize
    % Make an annulus for computing mean and variance
    %%%%% hardwired??
    ringo=0.48*n0;
    ringi=0.42*n0;
    msko=fuzzymask(n0,2,ringo,ringo/10);
    mski=fuzzymask(n0,2,ringi,ringi/10);
    annulus=msko-mski;
    
    % Compute mean and variance
    nann2=annulus(:)'*annulus(:);
    nann=sum(annulus(:));
    vars=zeros(nImgs,1);
    avgs=zeros(nImgs,1);
    for i=1:nImgs
        pix=allImgs(:,:,i).*annulus;
        avgs(i)=sum(pix(:))/nann;
        vars(i)=pix(:)'*pix(:)/nann2-avgs(i)^2;
    end;
    
    % Normalize variance to 1 in each image
    for i=1:nImgs
        imgs(:,:,i)=(allImgs(:,:,i)-avgs(i))/sqrt(vars(i));
    end;
    if readAltImgs % we'll use the same variance values.
        %         pix=allImgs(:,:,i).*annulus;
        %         avgs(i)=sum(pix(:))/nann;
        altImgs(:,:,i)=(allAltImgs(:,:,i)-avgs(i))/sqrt(vars(i));
    end;
end;

