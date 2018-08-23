% raLatentVars.m
% Look at roi.mat files and show the latent variables

    [roiFiles, pa]=uigetfile('*roi.mat','Select roi files','multiselect','on');
    if isnumeric(pa)
        return
    end;
    reconPath=AddSlash(pa);
    rootPath=ParsePath(ParsePath(pa));
    cd(reconPath);
    
    if ~iscell(roiFiles)
        roiFiles={roiFiles};
    end;
    
disp('Loading ri.mat');
ri=LoadStruct('ri.mat');

cd(rootPath)
disp('Loading stack');
si=LoadStruct([ri.stackPath ri.siName]);

refAngles=reGetAngleList(ri,false);  % pick beta, gamma angles

%%
nr=numel(roiFiles);
for i=1:1
    roiName=roiFiles{i};
    iTwin=double(roiName(4))-96 % a->1 , b->2
    iter=str2double(roiName(2:3))
    moiName=roiName;
    moiName(1,end-6)='m';
    cd(reconPath);
    disp(['Loading ' roiName]);
    roi=LoadStruct(roiName);
    disp(['Loading ' moiName]);
    moi=LoadStruct(moiName);

    
    disp('Loading images');
    ri.nGroups=1;
    cd(rootPath);
    [gSi,imgs]=reLoadStackGroup(ri,si,iTwin,1,0);

    disp('Making refs');
    refs=reMakeTemplates(moi.refVols,refAngles);  % refs(x,y,iRef,iVol)
%%
    imIndices=1:6;
    reShowLatentVars(imgs,refs,ri,moi,roi,iter,[0 1],imIndices)

end;    
    
return    
    
    disp('Making refs');
    refs=reMakeTemplates(moi.refVols,refAngles);  % refs(x,y,iRef,iVol)
    nim=size(roi.imgAmps,2);
    nVols=size(roi.pRefs,2);
    nRefs=size(roi.pRefs,1);
    nt=sqrt(size(roi.pTrans,1));
    prv=reshape(roi.pRefs,nVols*nRefs,nim); % lump together refs and volumes
    
%     Find the indices of best-match reference, translation and alpha
    [~,irv]=max(prv,[],1);
    [iRef,iVol]=ind2sub([nVols nRefs],irv);
    [~,it1]=max(roi.pTrans,[],1);
    [iTransx,iTransy]=ind2sub([nt nt],it1);
    [~,iAlpha]=max(roi.pAlphas,[],1);
 






nBin=3;
doCTF=1;
k=1e-7;

[names, pa]=uigetfile('*roi.mat','Select roi.mat files','multiselect','on');
if isnumeric(pa)
    return
end;
cd(pa);
if ~iscell(names)
    names={names};
end;

nFiles=numel(names);
fileIndex=1;
b=1;

while b ~= 'q'
    
    if fileIndex>nFiles
        fileIndex=nFiles;
        beep;
    elseif fileIndex<1
        fileIndex=1;
        beep;
    else
        roiName=names{fileIndex};
        
        % Load either roi or gRoi structure.
        disp(['loading ' roiName]);
        roi=loadStruct(roiName);
        moiName=roiName;
        
        
        % Show the class means
        [nx,ny,nim]=size(roi.classMeans);
        classMeansNorm=zeros(nx,ny,nim,'single');
        if doCTF
            for i=1:nim
                cNorm=ifftshift(roi.classNorms(:,:,i));
                if cNorm(:)'*cNorm(:)>k.^2 %% something there
                    classMeansNorm(:,:,i)=real(ifftn(fftn(roi.classMeans(:,:,i))./(k+cNorm)));
                end;
            end;
        else
            % Simple normalization by std
            cMeans=reshape(roi.classMeans,nx*ny,nim);
            std0=std(cMeans,1);
            sds=max(std0/1000,std(cMeans,1));
            classMeansNorm=reshape(cMeans./repmat(sds,nx*ny,1),nx,ny,nim);
        end;
        
        figure(1);
        ImagicDisplay3(BinImage(classMeansNorm,nBin));
        set(gcf,'name',roiName);
                
    end;  % if fileIndex in range
    
        
    b=0;
    while ~any(b=='npq')
        [partIndex,coords,b]=ImagicDisplay3;
        pause(0.1);
    end;
    
    switch b
        case 'n'
            fileIndex=fileIndex+1;
        case 'p'
            fileIndex=fileIndex-1;
    end;
end;
disp('Done.');


return


%%

% Look at the geometry of YClick vs rVesicle+-mbnOffset
si1=si;
np=numel(si1.miIndex);
rsos=false(np,1);
for i=1:np
    rsos(i)=si1.mi{si1.miIndex(i)}.particle.picks(si1.miParticle(i),7);
end;

sel=~rsos;
mo=si1.mbnOffset;
figure(3);
plot(si1.rVesicle(sel),si1.yClick(sel),'.',si1.rVesicle(sel),si1.rVesicle(sel)+mo,'-');
