% TestVesicleShifts
%
%  Make the modified mi file, and shift the micrographs.
basePath='/Volumes/TetraData/EMWork/Hideki/121210/TrackingTest/';
load([basePath 'Info/sq02_10000mi.mat']);
mi.basePath=basePath;
mi.imagePath='Micrograph/';
names=mi.imageFilenames;
nim=numel(names);
shifts=[ 10 0; 0 0; 0 10 ];
% Modify image files
for i=1:nim
    name=[mi.basePath mi.imagePath mi.imageFilenames{i}]
    m=ReadEMFile(name);
    msh=circshift(m,shifts(i,:));
    [pa nm]=fileparts(mi.imageFilenames{i});
    outName=[nm '.mrc'];
    mi.imageFilenames{i}=outName;
    fullOutName=[mi.basePath mi.imagePath mi.imageFilenames{i}]
    WriteMRC(msh,mi.pixA,fullOutName);
end;
save([basePath 'Info/sq02_10000miShifted.mat'],'mi');

%%

% Mark each particle position with a white dot.
% This doesn't work because I haven't transformed the images...

for i=1:nim
    name=[mi.basePath mi.imagePath mi.imageFilenames{i}]
    m=ReadEMFile(name);
    m=RemoveOutliers(m);
    mxm=max(m(:));
    np=size(mi.particle.picks,1);
    for j=1:np
        flag=mi.particle.picks(j,3);
        if flag>=16 && flag <=33
            coords=mi.particle.picks(j,1:2)+1;
            vind=mi.particle.picks(j,4);
            if vind>0
                shCoords=coords+[mi.vesicle.shiftX(vind,i) mi.vesicle.shiftY(vind,i)];
            else
                shCoords=coords;
            end;
            m(round(shCoords(1)),round(shCoords(2)))=mxm;  % mark a bright spot.
        end;
    end;
WriteMRC(m,mi.pixA,name);
end;
