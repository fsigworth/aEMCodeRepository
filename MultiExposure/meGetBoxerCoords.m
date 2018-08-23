% meGetBoxerCoords
% given paths, pick up mi files, see if there is a corresponding boxer
% file, and read it into the mi structure.  Also get the particle stack.

displayOn=1;
imIndex=2; % zero means merged.
boxSize=80;
imScale=2; % scale-up of picker coordinates; typical values are 1 or 2.

% baseDir='/Volumes/TetraData/EMWork/Hideki/110628/';
baseDir='/Volumes/TetraSSD/DSR from raid4/DSR wo carbon 2.5mgml-1 blot1sec/DDD/Images/';
cd(baseDir);
miDir='InfoCopy/';
miStr='mi.mat';
miStrLength=numel(miStr);
boxDir='box/';
boxStr='m.box';
if imIndex==0
    imageDir='Merged/';
    imageStr='m.mrc';
else
    imageDir='../../DDDImages/micrograph/';
end;

dmi=dir(miDir);
ndmi=numel(dmi);
dbox=dir(boxDir);
ndbox=numel(dbox);
stack=single(zeros(boxSize,boxSize,1));
nims=0;
j=0;
si=struct;

for i=1:ndmi
    miName=dmi(i).name;
    p=regexp(miName,miStr,'end');
    if numel(p)>0
        disp(miName);
        load([miDir miName]);
        boxName=[boxDir 'f' mi.baseFilename boxStr];
        if imIndex>0
            imageName=[imageDir mi.imageFilenames{imIndex}];
            tMat=(mi.mergeMatrix(:,:,imIndex));
        else
            imageName=[imageDir mi.baseFilename imageStr];
            tMat=eye(3);
        end;
        if FileExists(boxName) && FileExists(imageName)
            m0=ReadEMFile(imageName);
            if imIndex>0
                m0=RemoveOutliers(m0);
            end;
            ds=mi.imageSize(1)/size(m0,1);
            pixA=mi.pixA*ds;
            
            [imgs,pts]=ExtractBoxes(m0, boxName, boxSize,...
                imScale, displayOn, tMat);
            drawnow;
            nnew=size(imgs,3);
            startImageNo=nims+1;
            stack(:,:,nims+1:nims+nnew)=imgs;
            nims=nims+nnew;
            j=j+1;
            disp([imageName num2str(nnew)]);
            %             coordScale=mi.imageSize(1)/imgsize(1);  % scale up if we are picking from compressed micrographs
            mi=rmfield(mi,'particle');
            mi.particle.picks=single(zeros(nnew,8));
            mi.particle.picks(:,1:2)=pts';
            %             for j=1:nnew
            %                 mi.particle(j).x=pts(j,1);  % coords in original file
            %                 mi.particle(j).y=pts(j,2);
            %             end;
            % hold on;
            % plot(mi.particle.picks(:,1),mi.particle.picks(:,2),'r.');
            % hold off;
            save([miDir miName],'mi');
%             Get the CTF.  No flipping is performed if imIndex > 0.
% ***note that we now use the mergeMode variable
            if imIndex>0
                mergeMode=3;
            else
                mergeMode=1;
            end;
            si.ctfs(:,:,j)=meGetEffectiveCTF(mi,boxSize,ds,mergeMode);
            si.pixA=pixA;
            si.filename{j}=miName;
            si.mis{j}=mi;
            si.fileIndex(startImageNo:nims)=j;
            si.defocus(startImageNo:nims)=mi.ctf(max(imIndex,1)).defocus;
        end;
    end;
end;
%%
si
whos stack
outName=['Stack/DSR' num2str(imIndex)]
save([outName 'si.mat'],'si');
save([outName 'stack.mat'],'stack');

