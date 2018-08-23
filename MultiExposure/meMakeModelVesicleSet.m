function vs=meMakeModelVesicleSet(mi,mMerged,vindex,doCTF,extraFilter,nZeros,mode)
% function vs=meMakeModelVesicleSet(mi,mMerged,vindex,doCTF,extraFilter,nZeros,mode)
% function vs=meMakeModelVesicleSet(mi,n,vindex,doCTF,extraFilter,nZeros)
% Create a set of ctf-filtered vesicle fits, such that sum(vs,3) is an
% optimum fit to mMerged for subtraction.  (alternatively the 2nd argument
% can be just the size of the images in vs, in which case no ls fit is
% performed.  Only the first two arguments are required.
% If some of mi.weights are zero, the corresponding planes of vs will be
% zeros.
% extraFilter is a filter function applied to the output images, center
% zero frequency.
% mode is the merging mode; 1 is normal merging.

% test code
% load('/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/Info/sq02_10000mi.mat');
% mMerged=ReadMRC('/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/Merged/sq02_10000m.mrc');
% vindex=0;
% doCTF=1;
% doPW=0;
% nZeros=1;

% parameters to select 'good' vesicles to include in the LS fit.
maxSMult=2;
maxRTerm=5;


shiftMax=10;  % max pixels shifted

nv=numel(mi.vesicle.x);

if numel(mMerged)>2  % mMerged is an image
    n=size(mMerged);
    doLSFit=1;
else
    doLSFit=0;
    n=mMerged;       % mMerged is the size of the output
    if numel(n)<2
        n(2)=n(1);
    end;
end;

if nargin<3 || vindex(1)==0  % vindex=0 means compute all the vesicles.
    if ~isfield(mi.vesicle,'ok')
        mi.vesicle.ok=1:numel(mi.vesicle.x);
    end;
    vindex=find(all(mi.vesicle.ok,2));
else
    vindex = vindex(vindex<=nv);  % don't allow out-of-range indices
end;
if nargin<4
    doCTF=1;
end;
if nargin<5
    extraFilter=1;
end;
if nargin<6
    nZeros=1;
end;
if nargin<7
    mode=1;  % merging mode
end;

if numel(vindex)<1
    return
end;

nim=numel(mi.doses);
if numel(mi.weights)>nim
    mi.weights=mi.weights(1:nim);
end;
imSet=find(mi.weights~=0);

nima=numel(imSet);  % number of active images

ds=mi.imageSize(1)/n(1);   % downsample factor

vs=single(zeros([n nim]));  % default is a zero image.
vsGood=vs;
if doCTF
    freqs=RadiusNorm(n)/(mi.pixA*ds);
    %     [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs( freqs, mi.ctf, mi.doses, nZeros, mi.weights);
    [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs2( freqs, mi, nZeros, mode );
    
end;
mis(nim)=mi;  % force the size of the mi array
vs=single(zeros([n nim]));
oldTruncShiftx=inf;
oldTruncShifty=inf;


for i=imSet
    mis(i)=mi;
    %         Zero out all the weights except for the i'th one.
    mis(i).weights=0*mi.weights;
    mis(i).weights(i)=1;
    
    %         Shift the individual vesicle coordinates
    if isfield(mi.vesicle,'shiftX') && numel(mi.vesicle.shiftX)>0
        shiftx=mi.vesicle.shiftX(:,i);
        shifty=mi.vesicle.shiftY(:,i);
        truncShiftx=max(min(shiftx,shiftMax),-shiftMax);
        truncShifty=max(min(shifty,shiftMax),-shiftMax);
    else
        truncShiftx=0;
        truncShifty=0;
    end;
    
    mis(i).vesicle.x=mis(i).vesicle.x(:)+truncShiftx;
    mis(i).vesicle.y=mis(i).vesicle.y(:)+truncShifty;
    if any(truncShiftx~=oldTruncShiftx) || any(truncShifty~=oldTruncShifty)
        v0=meMakeModelVesicles(mis(i),n,vindex,0);  % no ctf
        if doLSFit  % we will pick only 'good' vesicles for ls fitting of amplitude.
            goodVes=all(mis(i).vesicle.ok(:,1:3),2)...
                & all(abs(mis(i).vesicle.r(:,2:end))<maxRTerm,2);  % exist, in range, refined
            distortedVes=any(abs(mis(i).vesicle.s(:,2:end))>maxSMult,2);
            nDistGood=sum(goodVes & distortedVes);
            if nDistGood>0
                disp([num2str(nDistGood) ' distorted vesicles excluded from vesicle fit.']);
                goodVes=goodVes & (~distortedVes);
            end;
            okVes=false(size(goodVes));
            okVes(vindex)=true;
            goodVes=goodVes & okVes;
            if sum(goodVes)<1
                warning('No good vesicles found for LS fitting');
                doLSFit=0;
            else % we will do least-squares only of good ves
                vGood=meMakeModelVesicles(mis(i),n,find(goodVes),0);
            end;
        end;
    end;
    if doCTF
        q=find(mi.doses,1);  % find the first nonzero dose
        if numel(q)<1, error('Doses are all zero'); end;
        dose1=mi.doses(q);
        % Filter the vesicles
        vFilter=ifftshift(extraFilter.*dctfs(:,:,i).*coeffs(:,:,i))...
            *mi.ctf(i).ampFactor*mi.doses(i)/dose1;
        vs(:,:,i)=real(ifftn(fftn(v0).*vFilter));
        if doLSFit % we will do least-squares only of good ves
            vsGood(:,:,i)=real(ifftn(fftn(vGood).*vFilter));
        end;
    else
        vs(:,:,i)=v0;
    end;
    %     subplot(2,2,i);
    %     imacs(vs(:,:,i));
    %     title(i);
    %     drawnow;
    oldTruncShiftx=truncShiftx;
    oldTruncShifty=truncShifty;
    
end;
%%
if doLSFit  % We have an image to fit to
    %     Do a least-squares fit of the good vesicles to the data image
    fs=reshape(vsGood(:,:,imSet),prod(n),nima);
    fs=[fs ones(prod(n),1)];
    as=LinLeastSquares(fs,mMerged(:));
    for j=1:nima
        i=imSet(j);
        vs(:,:,i)=as(j)*vs(:,:,i);
    end;
    disp(['Amp factors: ' num2str(as(:)')]);
end;

% subplot(2,2,4);
% imacs(mMerged-sum(vs,3));
% title(num2str(as'));
