% KvMakeMovieLinks.m
% make a directory full of links to the movie files.

rootDir= '~/tomo/161101/';
infoDirs={'KvLipo122_4a/' 'KvLipo122_4b/'};
movieDirs={'KvLipo122_4a_raw/' 'KvLipo122_4a_raw/'};
linkDir='Movies/';
infoDir='Info/';
    cd(rootDir);

% First, make a link for every movie known to an info file

for iInfo=1:numel(infoDirs)
    names=f2FindInfoFiles([infoDirs{iInfo} 'Info/']);
    for i=1:numel(names)
        mi=ReadMiFile(names{i});
        oldMvPath=[rootDir movieDirs{iInfo} mi.moviePath];
        disp([oldMvPath mi.movieFilename]);
        newMi=mi;
        newMi.moviePath=linkDir;
        newMi.movieFilename=[mi.baseFilename '.tif'];
        newMi.basePath=rootDir;
       linkName=[linkDir newMi.movieFilename];
        str=['ln -s ' oldMvPath mi.movieFilename ' ' linkName];
        disp(str)
        system(str);
        newMiName=[infoDir newMi.baseFilename '.mi.txt'];
        WriteMiFile(newMi,newMiName);
     end;
         
        
        
end;




