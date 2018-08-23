% meGetHidekiBoxes
% given paths, pick up mi files, see if there is a corresponding boxer
% file, and read it into the mi structure.  Also get the particle stack.

mode=2;  %Process the single-image (non-merged) files

displayOn=1;
figure(1); clf; SetGrayscale;
boxSize=64;
imScale=1; % scale of micrograph we're extracting from, relative to the one
% used for picking; typical values are 1 or 2.

switch mode
    case 1
        baseDir='/Volumes/TetraData/EMWork/Hideki/110628/';
        cd(baseDir);
        miDir='HS12slot4_GluR2-GFP/';
        miStr='mi.mat';
        miStrLength=numel(miStr);
        boxDir='box/';
        boxStr='merge.box';
        imageDir='Merged/';
        imageStr='m.mrc';
        imageDir=miDir;
        % 'MergedWhitened/';
        imageStr='mw.mrc';
        stackName='stack';
    case 2
        miDir='NoMerge/';
        miStr='mi.mat';
        miStrLength=numel(miStr);
        boxDir='box/';
        boxStr='merge.box';
        imageDir='NoMerge/';
        imageStr='m.mrc';
        imageDir=miDir;
        imageStr='mw.mrc';
        stackDir='';
        stackName='NoMergePWStack';
end;

dmi=dir(miDir);
ndmi=numel(dmi);
dbox=dir(boxDir);
ndbox=numel(dbox);
stack=single(zeros(boxSize,boxSize,1));
nims=0;
for i=1:ndmi
    miName=dmi(i).name;
    p=regexp(miName,miStr,'end');
    if numel(p)>0
        disp(miName);
        load([miDir miName]);
        boxName=[boxDir mi.baseFilename boxStr];
        imageName=[imageDir mi.baseFilename imageStr];
        if FileExists(boxName) && FileExists(imageName)
            [imgs,pts,imgsize]=ExtractBoxes(imageName, boxName, boxSize, imScale, displayOn);
            drawnow;
            nnew=size(imgs,3);
            stack(:,:,nims+1:nims+nnew)=imgs;
            nims=nims+nnew;
            disp([imageName num2str(nnew)]);
            coordScale=mi.nPix(1)/imgsize(1);  % scale up if we are picking from compressed micrographs
            for j=1:nnew
                mi.particle(j).x=pts(j,1)*coordScale;  % coords in original file
                mi.particle(j).y=pts(j,2)*coordScale;
            end;
            save([miDir miName],'mi');
        end;
    end;
end;

WriteImagic(stack,[stackDir stackName]);
