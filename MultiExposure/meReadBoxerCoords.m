% meReadBoxerCoords
% Reader Boxer db files and write the particle coordinates into mi files.
% Pick up mi files, see if there is a corresponding boxer
% file, and read it into the mi structure.  Show the image and overlaid
% boxes as a check.

displayOn=1;
boxSize=44;
writeOut=0;

replaceParticles=1;  % Remove any existing particles in mi files

baseDir='/Volumes/TetraData/EMWork/Hideki/121118/PC1PC2liposome_SUP_14k/';
miDir='Info/';
boxDir='Boxes/';

miStr='mi.mat';
boxStr='mf.box';
imageStr='m.mrc';

dmi=dir([baseDir miDir]);
ndmi=numel(dmi);
if ndmi<3
    error(['No mi files found in ' baseDir miDir]);
end;
dbox=dir([baseDir boxDir]);
ndbox=numel(dbox);
if ndbox<3
    error(['No box files found in ' baseDir boxDir]);
end;

for i=1:ndmi
    miName=dmi(i).name;
    p=regexp(miName,miStr,'end');
    if numel(p)>0
        disp(miName);
        load([baseDir miDir miName]);
        boxName=[baseDir boxDir mi.baseFilename boxStr];
        imageName=[baseDir mi.procPath mi.baseFilename imageStr];
        if FileExists(boxName) && FileExists(imageName)
            disp(['Reading box file ' boxName]);
            bfile=fopen(boxName);
            if bfile<1
                error(['invalid file name ' boxName]);
            end;
            bc=textscan(bfile,'%n %n %n %n %*[^\n]');
            fclose(bfile);
            pts=[bc{1} bc{2}]+[bc{3} bc{4}]/2+1;  % This is what EMAN boxer does.
            
            %%
            disp(['Reading image    ' imageName]);
            [img pixA]=ReadEMFile(imageName);
            [nx ny]=size(img);
            ds=mi.imageSize(1)/nx;
            [nim nc]=size(pts);
            
            ShowImageAndBoxes(GaussFilt(img,.1),pts,64);
            title(imageName,'interpreter','none');
            drawnow;
            pause
            %% put in the coordinates
            np0=numel(mi.particle.x);
            if np0>0
                disp([num2str(np0) ' particles already present']);
            end;
            if replaceParticles
                disp([num2str(nim) ' particles to write']);
                mi.particle.x=pts(:,1)*ds;
                mi.particle.y=pts(:,2)*ds;
                mi.particle.vesicle=0;  % null vesicle
                mi.particle.type=ones(np0,1);
            else
                p0=np0+1;
                disp([num2str(nim) ' particles to add']);
                mi.particle.x(p0:p0+nim-1)=pts(:,1)*ds;
                mi.particle.y(p0:p0+nim-1)=pts(:,2)*ds;
                mi.particle.vesicle(p0:p0+nim-1)=0;  % null vesicle
                mi.particle.type(p0:p0+nim-1)=1;
            end;
            if writeOut
                save([baseDir miDir miName],'mi');
            end;
        end;
    end;
end;

