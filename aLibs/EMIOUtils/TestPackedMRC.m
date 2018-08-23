% TestPackedMRC.m
% Test of the packed uint4 mode of WriteMRC and ReadMRC.
tic
disp('read original file');
q=ReadEMFile('Aug22_23.19.51.mrc');

WriteMRC(q,1,'xtestFull.mrc',0);
% q(end,:,:)=[];  % make odd line length for test.
whos q

toc
disp('write');
tic
WriteMRC(q,1,'xtest2.mrc',32);
toc
disp('read');
tic
r=ReadMRC('xtest2.mrc');
toc
p=r-min(q,15);
min(p(:))
max(p(:))



return


%%
q1=min(15,q);


figure(1);
hist(single(q1(:)),1000);
figure(2);
hist(single(r(:)),1000);
