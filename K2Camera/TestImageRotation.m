
% Read from my mrc
cd('/Users/fred/EMWork/Hideki/140516/KvBetaLiposome_pH8_washedpH8/movie_frames/sq03_1')
mrf=RemoveOutliers(sum(ReadMovie('May16_18.24.57.mrc'),3));

% Read DM's TIFF
[m2, s2]=ReadMovie('May16_14.43.59.tif');
m2f=RemoveOutliers(sum(m2,3));

ref2=RemoveOutliers(ReadEMFile('CountRef_May16_11.27.37.dm4'));

cd('/Users/fred/EMWork/Hideki/140105/Box_AMPAR_Liposome_26_10mMGlu_slot4/movie_frames/sq02_0')
[m1,s1]=ReadMovie('Jan05_12.31.03.tif');
m1f=RemoveOutliers(sum(m1,3));

ref1=RemoveOutliers(ReadEMFile('CountRef_Jan05_12.31.03.dm4'));

%%
figure(1); SetGrayscale;
subplot(2,2,1);
imac(imscale(GaussFilt(m1f.*rot90(ref1),.02),256,[.03 0]));
title('My TIFF or MRC');
drawnow;

subplot(2,2,2);
imacs(GaussFilt(m2f.*rot90(ref2),.02));
title('TIFF from SerialEM');

subplot(2,2,3);
imacs(rot90(GaussFilt(ref1,.02)));
title('Ref rotated 90 degrees');

subplot(2,2,4);
imacs(rot90(GaussFilt(ref2,.02)));


