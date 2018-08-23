% LoadDSRStack

name='/Volumes/TetraData/EMWork/Hideki/DSR_Dataset/start';

m=ReadImagic(name);
%%
mf=GaussFilt(m,.3,1);
%%
[sp mfn]=spNormalizeAndCenter(mf,5.9,.52);

rawStack=uint8(40*Downsample(mfn,32,1)+128);

imats(rawStack,2);

save stackDSR32 rawStack

