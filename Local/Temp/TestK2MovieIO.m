% TestK2MovieIO

cd('/Volumes/WDTera2/EMWork/Hideki/140618/KvBetaLiposome_x3_long_movies/movies_long')
tic
disp('Read compressed TIFF');
m=ReadMovie('Jun18_15.23.05.tif');
toc

sz=size(m);
disp(['Stack size: ' num2str(sz)]);
disp(' ');

tic
disp('Convert to float');
fm=single(m);  % convert to float
toc
disp(' ');

tic
disp('Convert back to byte');
m2=uint8(fm);
toc;
disp(' ');

tic
disp('Write 8-bit MRC to USB disk')
WriteMRC(m,1,'byteTest.mrc',0);
toc
disp(' ');

tic
disp('Read 8-bit MRC from USB disk');
mr=ReadMRC('byteTest.mrc');
toc
disp(' ');

tic
disp('Write 32-bit MRC to USB disk');
WriteMRC(fm,1,'floatTest.mrc',2);
toc
disp(' ');

%%
cd ~/EMWork/Hideki/XX
tic
disp('Write 8-bit MRC to SSD');
WriteMRC(m,1,'byteTest.mrc',0);
toc
disp(' ');

tic
disp('Write 32-bit MRC to SSD');
WriteMRC(fm,1,'floatTest.mrc',2);
toc
disp(' ');
