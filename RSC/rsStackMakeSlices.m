function sis=rsStackMakeSlices(si,imgsPerSlice,nTwins,activeFlags)
% function sis=rsStackMakeSlices(si,imgsPerSlice,nTwins,activeFlags)
% Takes a StackInfo structure and creates an array of such structures sis,
% ns x nTwins in size, where ns is the number of slices.
% The output sis(i,j) has the additional field imgRange.  To make the
% appropriate local image stack, use
%   ir=sis(i,j).imgRange;
%   localImgs=imgs(:,:,ir(1):ir(2):ir(3));
% activeFlags is optional array of booleans, one per image.

nim=numel(si.miIndex);
if nargin<4
    activeFlags=true(nim,1);
end;
if imgsPerSlice<2
    imgsPerSlice=1;
end;

activeInds=find(activeFlags);
nim=numel(activeInds);
nSlices=ceil(nim/nTwins/imgsPerSlice);


sis=struct;
sis(nSlices,nTwins).pixA=0;  % allocate the array
si1=struct;
si1.pixA=si.pixA;

for i=1:nSlices
    for j=1:nTwins
        istart=(i-1)*imgsPerSlice+j;
        iend=min(i*imgsPerSlice+j-nTwins,nim);
        inds=(istart:nTwins:iend)';
        imgInds=activeInds(inds);
        
        activeMis=single(si.miIndex(imgInds));
        h=hist(activeMis,1:numel(si.mi));  % find out which mi's are used
        miFlags=(h>0)';
        miPtrs=cumsum(miFlags);  % lookup global to local pointers
        miInds=find(miFlags);    % local ponters to global 
        
        si1.miIndex=miPtrs(si.miIndex(imgInds));
        si1.miParticle=si.miParticle(imgInds);
        si1.alpha0=si.alpha0(imgInds);
        si1.rVesicle=si.rVesicle(imgInds);
        si1.sVesicle=si.sVesicle(imgInds);
        si1.mi=si.mi(miInds);
        si1.ctfs=si.ctfs(:,:,miInds);
        si1.imgRange=[istart nTwins iend];
        if i==1 && j==1
            sis=si1;
        else
            sis(i,j)=si1;
        end;
    end;
end;
