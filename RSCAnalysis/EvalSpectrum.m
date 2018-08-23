cd('/Users/fred/EMWork/Hideki/151216/Info')
% name=uigetfile
mi=ReadMiFile('g5_0010_Dec16_16.08.36mi.txt');
% mi=ReadMiFile('g5_0337_Dec17_00.31.44mi.txt');
m=meReadMergedImage(mi);
[s,freqs]=RadialSpectrum(m,mi.ctf(1),mi.pixA*2);
mi.ctf(1).B=10;
c2=meGetEffectiveCTF(mi,size(m)).^2;
c1=RadialCTF(c2,mi.ctf(1),mi.pixA*2);
plot(freqs,[GaussFilt(min(.3,s),.1) c1*.07+.03]);
title(['151215 ' mi.baseFilename ' B=10'],'interpreter','none');
ylabel('Spectral density of normalized image');
xlabel('Frequency, Å^{-1}');
