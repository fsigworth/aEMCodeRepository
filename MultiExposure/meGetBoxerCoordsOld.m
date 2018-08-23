% meReadBoxerCoords
% given paths, pick up mi files, see if there is a corresponding boxer
% file, and read it into the mi structure.  Also get the particle stack.

displayOn=1;
boxSize=44;
imScale=1; % scale of micrograph we're extracting from, relative to the one
% used for picking; typical values are 1 or 2.

baseDir='/Volumes/TetraData/EMWork/Hideki/121118/PC1PC2liposome_SUP_14kcd/';
miDir='Info/';
miStr='mi.mat';
boxDir='Boxes/';
boxStr='merge.box';
boxStr='mf.box';
% imageDir='MergedWhitened/';
miStrLength=numel(miStr);
imageStr='m.mrc';

dmi=dir(miDir);
ndmi=numel(dmi);
dbox=dir(boxDir);
ndbox=numel(dbox);
nims=0;
for i=1:ndmi
    miName=dmi(i).name;
    p=regexp(miName,miStr,'end');
    if numel(p)>0
        disp(miName);
        load([miDir miName]);
        boxName=[boxDir mi.baseFilename boxStr];
        imageName=[mi.procPath mi.baseFilename imageStr];
        if FileExists(boxName) && FileExists(imageName)
            imageName
            disp(['Reading box file ' boxName]);
            bfile=fopen(boxName);
            if bfile<1
                error(['invalid file name ' boxName]);
            end;
            bc=textscan(bfile,'%n %n %n %n %*[^\n]');
            fclose(bfile);
            %%
            disp(['Reading image    ' imageName]);
            [img pixA]=ReadEMFile(imageName);
            [nx ny]=size(img);
            % img=rot90(single(img-mean(img(:))),3); %!!!!!!!!
            % Compute the centers of the boxes
            % pts=[nx-bc{2}-bc{4}/2+1 bc{1}+bc{3}/2+1];
            pts=[bc{1} bc{2}]+[bc{3} bc{4}]/2+1;  % This is what EMAN boxer does.
            % pts=[nx-bc{2} bc{1}];
            % pts=[bc{1} nx-bc{2}];
            % pts=[bc{2} nx-bc{1}];
            [nim nc]=size(pts);
            ShowImageAndBoxes(GaussFilt(rot90(img,0),.1),pts,64);
            title(imageName,'interpreter','none');
            drawnow;
            %             [imgs,pts,imgsize]=ExtractBoxes(imageName, boxName, boxSize, imScale, displayOn);
            %             drawnow;
            %             nnew=size(imgs,3);
            %             stack(:,:,nims+1:nims+nnew)=imgs;
            %             nims=nims+nnew;
            %             disp([imageName num2str(nnew)]);
            %             coordScale=mi.nPix(1)/imgsize(1);  % scale up if we are picking from compressed micrographs
            %             for j=1:nnew
            %                 mi.particle(j).x=pts(j,1)*coordScale;  % coords in original file
            %                 mi.particle(j).y=pts(j,2)*coordScale;
            %             end;
            %             save([miDir miName],'mi');
        end;
    end;
end;

