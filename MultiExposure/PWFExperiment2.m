% PWFExperiment2

imagePath='/Volumes/TetraData/EMWork/Scott/';
 imageName='11oct20b_20170729_04_4096x4096_bright_0.mrc';
imageName='carbon.mrc';
% imageName='floodbeam.mrc';
% imageName='10sep19a_a_00006gr_00022sq_v01_00002hl_v01_00003en.mrc';
name=[imagePath imageName];
outname=[name 'o'];

CompressMRCFile(name,outname,2);

