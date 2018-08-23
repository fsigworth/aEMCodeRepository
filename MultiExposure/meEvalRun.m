% meEvalRun
ds=2;
i1=1;
i2=2;
m1=AffineTransform(BinImage(m(:,:,i1),ds),(mi.mergeMatrix(:,:,i1)));
m2=AffineTransform(BinImage(m(:,:,i2),ds),(mi.mergeMatrix(:,:,i2)));
pixA=mi.pixA*ds;
ndis=128;
fig=5;
CPars=mi.ctf(i1);
CPars(2)=mi.ctf(i2);
npanels=8;
meEvalLocalCorrelations2(m1,m2,pixA,CPars,npanels,ndis,fig);
