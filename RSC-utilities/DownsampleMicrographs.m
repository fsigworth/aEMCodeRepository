% DownsampleMicrographs
ds=4;
infoDir='Info/';
load([infoDir 'allMis.mat']);
nmi=numel(allMis);
outPath=('Micrograph_ds/');
CheckAndMakeDir('Micrograph_ds',1);
for i=1:nmi
    mi=allMis{i};
    name=[mi.imagePath mi.imageFilenames{1}];
    [pa,nm,ex]=fileparts(name);
    outName=[outPath nm 's' ex];
    m=ReadMRC(name);
    n=size(m);
    ms=Downsample(m,n/ds);
    disp(outName);
    WriteMRC(ms,mi.pixA*ds,outName);
end;    
