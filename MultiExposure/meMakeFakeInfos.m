% meMakeFakeInfos
% put up a file selector to get the base path, and then a file selector to
% get a representative image.
% Create an Info directory (in the base directory) and write an info file
% for each image.
% This allows for example jpeg merged images to be analyzed.

defoci=[1 10];
doses=[10 10];
pixA=5.8;
kV=200;

mi=meCreateMicrographInfoStruct11;
mi.pixA=pixA;
mi.keV=kV;
mi.doses=doses;
mi.ctf=meSetCTFitPars(defoci,pixA,kV);
for i=1:numel(defoci)
    mi.ctf(i).defocus=defoci(i);
    mi.ctf(i).deltadef=0;
end;

pa=uigetdir(pwd,'Select the base directory');
cd(pa)
mi.basePath=AddSlash(pwd);
[name imagePath]=uigetfile('*','Select an image file');
if numel(strfind(imagePath,mi.basePath))>0 % basePath is a substring
    mi.imagePath=AddSlash(imagePath(numel(mi.basePath)+1:numel(imagePath)));
    mi.procPath=mi.imagePath;
else
    error('Base path is not a substring of info path');
end;
%%
m=ReadEMFile([mi.procPath name]);
mi.imageSize=size(m);

mi.infoPath='Info/';
if ~DirectoryExists(mi.infoPath)
    mkdir(mi.infoPath);
end;
rexp='(?i).+\.(mrc|dm3|tif|jpg)';
names=FindFilenames(imagePath,rexp,0);
char(names)
nim=numel(names);
for i=1:nim
    [pa base ext]=fileparts(names{i});
    base=base(1:numel(base)-1);
    mi.baseFilename=base;
    mi.imageFilenames=names(i);
    save([mi.infoPath mi.baseFilename 'mi.mat'],'mi');
end;

