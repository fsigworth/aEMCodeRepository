% rtTiffMoviesToMRC

sourceDir='movie_frames/';
targetDir='mrc_movie_frames/';

CheckAndMakeDir(targetDir,1);

names=rtGetFilenames(sourceDir,'.tif',0);
disp([num2str(numel(names)) ' files.']);
for i=1:numel(names)
    disp(names{i});
    [mv,s]=ReadMovie([sourceDir names{i}]);
    [pa,nm,ex]=fileparts(names{i});
    WriteMRC(mv,s.pixA,[targetDir nm '.mrc'],0);
end;
