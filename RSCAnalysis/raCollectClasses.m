% raCollectClasses.m
% From roi structures, assign classes to the particles.
rootPath=uigetdir('','Select the experiment directory');
cd(rootPath);

    [roiFiles, roiPa]=uigetfile('*roi.mat','Select roi files','multiselect','on');
    if isnumeric(roiPa)
        return
    end;
    roiPa=AddSlash(roiPa);
    
    if ~iscell(roiFiles)
        roiFiles={roiFiles};
    end;
    
disp('Loading ri.mat');
ri=LoadStruct([roiPa 'ri.mat']);

disp('Loading si file');
si=LoadStruct([ri.stackPath ri.siName]);

refAngles=reGetAngleList(ri,false);  % pick beta, gamma angles
%%
nr=numel(roiFiles);
for i=1:nr
    roiName=roiFiles{i};
    disp(['Loading ' roiName]);
    roi=LoadStruct([roiPa roiName]);
    
    n=size(roi.classMeans,1);
    disp(['Scaling to n = ' num2str(n)]);
    ri=reMakeRunInfoScaled(ri,n);
    
    moiName=roiName;
    moiName(1,end-6)='m';
    disp(['Loading ' moiName]);
    moi=LoadStruct([roiPa moiName]);

    disp('Making refs');
    refs=reMakeTemplates(moi.refVols,refAngles);  % refs(x,y,iRef,iVol)
    
    iTwin=double(roiName(4)-'a')+1;
    if abs(iTwin-1.5)>.5
        error(['Couldn''t deduce iTwin from name: ' roiName]);
    end;
    iTwin
    ri.nGroups=1;
    disp('loading images');
    [gSi,gImgs]=reLoadStackGroup(ri,si,iTwin,1,0);

    disp('Collecting image statistics.');
    nim=size(roi.imgAmps,2);
    nVols=size(roi.pRefs,2);
    nRefs=size(roi.pRefs,1);
    nt=sqrt(size(roi.pTrans,1));
    prv=reshape(roi.pRefs,nVols*nRefs,nim); % lump together refs and volumes
    
%     Find the indices of best-match reference, translation and alpha
    [~,irv]=max(prv,[],1);
    [iRef,iVol]=ind2sub([nRefs nVols],irv);
    [~,it1]=max(roi.pTrans,[],1);
    [iTransx,iTransy]=ind2sub([nt nt],it1);
    [~,iAlpha]=max(roi.pAlphas,[],1);
    
%     format for writing
    gPi.iAlpha=int16(iAlpha);
    gPi.iRef=int16(iRef);
    gPi.iTransX=single(iTransx);
    gPi.iTransY=single(iTransy);
    gPi.refAlphas=ri.alphasI;
    gPi.refAngles=ri.angles;
    
    outName=[roiPa roiName(1:4)];
    disp(['Writing output files ' outName]);
    save([outName 'gImgs.mat'],'gImgs');
    save([outName 'gSi.mat'],'gSi');
    save([outName 'gPi.mat'],'gPi');
    
%     Do the translation and rotation of images
    disp('Rotating images');
    cImgs=rsRotateImage(gImgs,-gPi.refAlphas(gPi.iAlpha));
    nim=size(cImgs,3);
    for j=1:nim
        cImgs(:,:,j)=circshift(cImgs(:,:,j),8-[gPi.iTransX(j) gPi.iTransY(j)]);
    end;
%     Saving aligned images
    save([outName 'cImgs.mat'],'cImgs');
end;

