% GetStackCTFs
% Get the effective CTF from each info file in a directory, and write a
% corresponding .mat file in a second directory.
% 
outputDir='ctfs/';
n=64;  % output size
ds=1;  % downsampling factor of images

infoDir=AddSlash(uigetdir(pwd,'Get a directory of info files'));
baseDir=ParsePath(infoDir);
cd(baseDir);
%%
infos=meLoadInfoFiles(infoDir);
ni=numel(infos);
%%
for i=1:ni
    mi=infos{i};
    ctf=meGetEffectiveCTF(mi,n,ds);
    ctName=[outputDir mi.baseFilename 'mctf.mat'];
    save(ctName,'ctf');
end;
    
