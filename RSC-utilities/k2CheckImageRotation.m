% K2CheckImageRotation
movieDir='movie_frames/sq10/ref/';
% movieDir='movie_frames/sq10/without_rotation_flip_MRC/';
refFile='movie_frames/sq07/CountRef_Jun25_11.50.49.dm4';
mRef=ReadEMFile(refFile);
rexp='.*\.tif';
% rexp='.*\.mrc';

figure(1);
SetGrayscale;
subplot(1,2,2);
imac(imscale(-mRef,256,.01));
mr=BinImage(mRef,8);
drawnow;

names={};
j=0;
d=dir(movieDir);
msum=[];
for i=1:numel(d)
    q=regexp(d(i).name,rexp);
    if ~d(i).isdir && numel(q)>0
        j=j+1;
        names{j}=d(i).name;
    end;
end;
names'

numFiles=j;
%%
imgs={};
for j=1:numFiles
    [m,pixA]=ReadMovie([movieDir names{j}]);
    if j==1
        msum=sum(single(m),3);
    else
        msum=msum+sum(single(m),3);
    end;
    subplot(1,2,1);
     mb=BinImage(msum,8);
%      mb=flip(BinImage(msum,8),2);
     imacs(mb);
%     imacs(mb.*mr);
    title(names{j},'interpreter','none');
    drawnow;
end;


