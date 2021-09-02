% HRPickingTest.m
% Locate particles on micrographs
skipStar=1;
if ~skipStar
    % cd ('~/EMWork/20210224_YW_sel/');
    dataStarName='Refine3D/job110/run_data.star';
    
    disp(['Reading ' dataStarName '...']);
    [nms,dat]=ReadStarFile(dataStarName);
    op=dat{1};
    d=dat{2};
    disp('Decoding image names...');
    [d.imgNos,d.imgFiles]=rlDecodeImageName(d.rlnImageName);
    %%
    np=numel(d.rlnImageName);
    [stackUNames,stackUPtrs,stackRPtrs]=unique(d.imgFiles);
    [micUNames,micUPtrs,micRPtrs]=unique(d.rlnMicrographName);
    nU=numel(micUNames);
end;

%% Gather a stack of particles
stackInd=1;
stackName=stackUNames{stackInd};
disp(stackName);
[stack,sts]=ReadMRC(stackName);
stackInds=find(stackRPtrs==stackInd);
% stackInds=1:15; %%%%
n=size(stack,1);
np=numel(stackInds);

%% Get our 3D ref
[tmRef,s]=ReadMRC('HRPicking/tmMap.mrc');
% micRef=ReadMRC('HRPicking/micMap.mrc');
fc=.25;
% tmRef1=SharpFilt(Crop(tmRef,n),fc/s.pixA);
% tmRef1=GaussFilt(Crop(tmRef,n),.2);
tmRef1=Crop(tmRef,n);

%% Get a ctf
ct=rlStarLinesToCtf(nms,dat,1);
[projs,angs,shifts]=rlMakeRelionProjections(tmRef1,d,stackInds,ct.pixA,10);
angs
% projs=squeeze(sum(tmRef1,1));
%%
ct.B=50;
ctProjs=zeros(n,n,np,'single');
c=ifftshift(CTF(n,ct));
for i=1:np
%      ctProjs(:,:,i)= real(ifftn(fftn(projs(:,:,i)).*c));
     ctProjs(:,:,i)= real(ifftn(fftn(projs(:,:,i)).*abs(c)));
%     ctProjs(:,:,i)=-GaussFilt(projs(:,:,i),.1) ;
end;

for i=1:np 
    j=stackInds(i);
subplot(222);
imags(-ctProjs(:,:,i));
subplot(221);
imgFlip=real(ifftn(fftn(stack(:,:,j)).*(sign(c))));
% imags(GaussFilt(stack(:,:,j),.1));
imags(GaussFilt(imgFlip,.1));
title([num2str(i) ' : ' num2str(angs(i,:),3)]);
subplot(223);
imags(fftshift(real(ifftn(fftn(stack(:,:,j)).*conj(fftn(-ctProjs(:,:,i)))))));
subplot(224);
imags(GaussFilt(projs(:,:,j),.1));
pause;
end;
return








%%
for i=1:100;
    if ~exist(micUNames{i},'file')
        continue;
    end;
    micPtrs=find(micRPtrs==i);
    np=numel(micPtrs);
    %     if numel(micPtrs)>0
    [mic,s]=ReadMRC(micUNames{i});
    fMic=GaussFilt(mic,.05);
    [bX,bY,tX,tY]=MakeBoxDrawingVectors( ...
        [d.rlnCoordinateX(micPtrs) d.rlnCoordinateY(micPtrs)], ...
        op.rlnImageSize(1)/2,0.8);
    tStrings=cell(np,1);
    for j=1:np
        tStrings{j}=num2str(micPtrs(j));
    end;
    imags(fMic);
    hold on;
    plot(bX,bY,'color',[1 1 0]);
    text(tX,tY,tStrings,'color',[1 1 0], ...
        'HorizontalAlignment','left',  'VerticalAlignment','top');
    hold off;
    title(i);
    drawnow;
end;

