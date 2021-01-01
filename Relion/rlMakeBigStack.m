% rlMakeBigStack.m

% Read the particles.star file from an Extract job. Read each particle
% image and assemble it into a stack. Write out the stack and a
% particle_stack.star file.

inDir='Extract/job066/';  % get the non-normalized stack
inStar='particles.star';
outStackName='BigStackSub2.mrcs';
outStar='particles_BigStackSub2.star';

disp(['Reading ' inDir inStar]);
[nms,dat]=ReadStarFile([inDir inStar]);
pixA=dat{1}.rlnImagePixelSize;
d=dat{2};
%%
nl=numel(d.rlnImageName);
newPartNos=zeros(nl,1);
newImgName=cell(nl,1);
newMicrographName=cell(nl,1);
currentFile='';
currentGroup='';
nStack=0;
for i=1:nl
    [partNo,fileName]=rlDecodeImageName(d.rlnImageName{i});
    if ~strcmp(fileName,currentFile) % Read a new file
        currentFile=fileName;
        imgs=ReadMRC(fileName); % imgs are type single
        nImgs=size(imgs,3);
        disp([fileName ' ' num2str(nImgs)]);
        if i==1
            n=size(imgs,1);
            bigStack=zeros(n,n,nl,'single');
        end;
            bigStack(:,:,nStack+1:nStack+nImgs)=imgs;
        offset=nStack;
        nStack=nStack+nImgs;
    end;
    if ~strcmp(d.rlnGroupName{i},currentGroup)
        currentGroup=d.rlnGroupName{i};
        micNameCell=d.rlnMicrographName(i);
    end;
    newPartNo=partNo+offset;
    newImgName{i}=[num2str(newPartNo,'%6d') '@' outStackName]; 
    newPartNos(i)=newPartNo;
    newMicrographName(i)=micNameCell;
end;      
if nStack<nl
    disp([whos bigStack
        nl
        'Removing excess bigStack planes: ' num2str([nStack nl])]);
    bigStack(:,:,nStack+1:nl)=[];
end;
%%
% Trim distributions of mean and sd by setting active flags
trimThresh=.001;
bigStack=reshape(bigStack,n^2,nStack);
means=mean(bigStack);
refMean=median(means);
upperMean=Percentile(means,1-trimThresh);
okMeans=means<upperMean & means>2*refMean-upperMean;
disp([num2str(sum(~okMeans)) ' out of range means, of ' num2str(nStack)]);

sds=std(bigStack);
refSD=median(sds);
lowerSD=Percentile(sds,trimThresh);
okSDs=sds>lowerSD & sds<2*refSD-lowerSD;
disp([num2str(sum(~okSDs)) ' out of range sds']);
activeFlags=okMeans & okSDs;

figure(2);
clf;
subplot(221);
nBins=ceil(sqrt(nStack));
hist(means(activeFlags),nBins);
title('Particle image means');
ylabel('Before normalization');

subplot(222);
hist(sds(activeFlags),nBins);
title('Standard Deviations');
drawnow;

%%
bigStack=reshape(bigStack,n,n,nStack); % restore the original size.

disp('Normalizing the stack');
nBlock=5000; % Do this in blocks to save memory.
istart=1;
while istart<nStack
    iend=min(istart+nBlock-1,nStack);
    bigStack(:,:,istart:iend)=-NormalizeImages(bigStack(:,:,istart:iend),1,0);
    istart=istart+nBlock;
    fprintf('.');
end;
fprintf('\n');
% Look at mean and sd after normalization
bigStack=reshape(bigStack,n^2,nStack);
means=-mean(bigStack);
refMean=median(means);
upperMean=Percentile(means,1-trimThresh);
okMeans2=means<upperMean & means>2*refMean-upperMean;

sds=std(bigStack);
refSD=median(sds);
lowerSD=Percentile(sds,trimThresh);
okSDs2=sds>lowerSD & sds<2*refSD-lowerSD;
ok2=okMeans2 & okSDs2;
disp([num2str(sum(~ok2)) ' additional out of range']);
activeFlags=activeFlags & ok2;

bigStack=reshape(bigStack,n,n,nStack); % restore the original size.

subplot(223);
hist(means(activeFlags),nBins);
ylabel('After normalization');
subplot(224);
hist(sds(activeFlags),nBins);
drawnow;

%%
%
disp(['Writing ' outStackName]);
WriteMRC(bigStack,pixA,outStackName);
%%

disp(['Writing ' outStar]);
dat{2}.rlnImageName=newImgName;
dat{2}.rlnMicrographName=newMicrographName;
WriteStarFile(nms,dat,outStar,'',{[] activeFlags});
disp('done.');



function minX=FindThreshVal(x,fraction);
nx=numel(x);
[xSort,xInds]=sort(x);
minInd=find(xSort>nx*thresh);
boolInds=false(nx,1);
boolInds(minInd:end)=true;
aboveThresh=boolInds(xInds);
end