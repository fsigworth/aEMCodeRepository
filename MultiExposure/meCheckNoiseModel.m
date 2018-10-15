% meCheckNoiseModel

cd '/Users/fred/EMWork/Hideki/160909p/KvLipo121_2w11v3m3'
mi=ReadMiFile('Info/sq02_1_0001_Sep09_17.23.21mi.txt');
mv=ReadMRC(['Merged/' mi.baseFilename 'mv.mrc']);
ds=4;
n=mi.imageSize(1)/ds;
mvs=Downsample(mv,n);
sp1=RadialPowerSpectrum(mvs)*mi.pixA^2*ds^2;
%%
f=(0:n/2-1)'/(mi.pixA*n);
[spec,shot]=meEvalNoiseModel(f,mi);
c2d=meGetEffectiveCTF(mi,n);
c=sectr(c2d);
semilogy(f,[spec.*(c.^2)+shot shot sp1])
