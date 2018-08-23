function TestParallelTiffReading

cd('/Users/fred/EMWork/Hideki/160909/KvLip121_1/movie_frames/sq04_1')
fname='Sep09_14.57.51.tif';
nim=30;
m=zeros(15,3710,3838,2,'uint8');
tic
parfor i=1:15
    m(i,:,:,:)=ReadTiffStack(fname,i,2);
end
size(m);
toc