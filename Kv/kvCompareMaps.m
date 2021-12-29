% kvCompareMaps.m
% Compare Kv structures in the 20211122 dataset

cd('/Users/fred/EMWork/Yangyu/20211122_farnam');

m0=ReadMRC('PostProcess/job088/postprocess.mrc');
m0r=rsRotateImage(m0,15);
m0x=Downsample(m0r,256,1);
m0x=Crop(m0x,200,1);

m1=ReadMRC('PostProcess/job057/postprocess.mrc');
m1d=Downsample(m1,128);
m1f=GaussFilt(m1d,.15);
m1x=Downsample(m1f,256,1);
m1x=Crop(m1x,200,1);

%%

figure(1);
imagsar(m1x(:,:,108:-1:1));



figure(2);
imagsar(m0x(:,:,108:-1:1));
