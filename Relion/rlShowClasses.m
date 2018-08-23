% rlShowClasses
% Put up a file selector for *si.mat,
disp('Getting a model.star file.');
[stName, stPath]=uigetfile('*model.star','Select star file');
if isnumeric(stPath)  % user has clicked Cancel
    return
else
    cd(stPath);
end;
[stBlocks,stData,ok]=ReadStarFile(stName);
% The classes are the second block in the file
gen=stData{1};
n=gen.rlnOriginalImageSize;
pixA=gen.rlnPixelSize;
ncls=gen.rlnNrClasses;
imgs=zeros(n,n,ncls,'single');
p=zeros(ncls,1);

cls=stData{2};
%%
doSort=1;

for i=1:ncls % loop through each of the images
    [ind,name]=rlDecodeImageName(cls.rlnReferenceImage{i});
[pa nm ex]=fileparts(name);
xName=[nm ex];
    imgs(:,:,i)=ReadMRC(xName,ind,1);
    p(i)=cls.rlnClassDistribution(i);
end;

if doSort
[pSort,iSort]=sort(p);
else
    pSort=p;
    iSort=1:ncls;
end;
nc=10;
nr=ceil(ncls/nc);
for i=1:ncls
    j=iSort(i);
    mysubplot(nr,nc,i);
    imags(fliplr(imgs(:,:,j)));
    axis off;
    title([num2str(j) '    ' num2str(100*p(j),3)]);
end;

