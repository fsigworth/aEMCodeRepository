% nSumF2SFrames

frames=1:5;
mvPath='../Sums13_6833/';
%load([mvPath 'allNames.mat']);

nmi=numel(allNames);
disp([num2str(nmi) ' files']);
for i=1:nmi
    mi=ReadMiFile(allNames{i});
    mvName=[mvPath mi.imagePath mi.baseFilename 'ali_Stk.mrc'];
 %   disp(mvName);;
    mv=ReadMRC(mvName);
    img=sum(single(mv(:,:,frames)),3);
    imgName=[mi.imagePath mi.imageFilenames{1}];
    disp([num2str(i) ' ' imgName]);
    WriteMRC(img,mi.pixA,imgName);
end;
