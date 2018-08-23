% EvalSpectrum2.m
% Do quantitative analysis of spectrum of a movie.
% In version 2 we re-align the movie with no damage compensation, then compute the radial
% spectrum.

doAlignment=0;

cd('/Users/fred/EMWork/Hideki/151216')
names{1}='Info/g5_0017_Dec16_16.25.41mi.txt';
% names{2}='Info/g5_0021_Dec16_16.34.58mi.txt';  no good Krios exp
% names{2}='Info/g5_0023_Dec16_16.38.18mi.txt';  ditto
names{2}='Info/g5_0024_Dec16_16.39.25mi.txt';
names{3}='Info/g5_0337_Dec17_00.31.44mi.txt';

if doAlignment
    tpars=struct;
    tpars.nameSegment='i';  % micrograph will be 'ali'
    tpars.forceFitting=0;   % use the existing shift information
    tpars.doDamageComp=0;
    tpars.doSaveFigures=0;
    tpars.writeZTiff=1;
    % Special for OHUS Krios data:
    tpars.oldK2=1;
    tpars.finalRot90=3;
    for i=1:3
        k2DriftTracker(names(i),tpars);
    end;
end;
return

%%
figure(2);
i=2
%     mpars=struct;
%     mpars.overwrite=1;
%     mpars.writeZTiff=1;
%     mpars.mergeMode=4;  % no damage compensation
%     mpars.ctfOptions.defoci=1;
%     MergeImages(names(i),mpars);

    mi=ReadMiFile(names{i});
%     imgName=[mi.imagePath mi.imageFilenames{1}];
%     [imgName,ok]=CheckForImageOrZTiff(imgName);
% if ok
%         m=ReadEMFile(imgName);
%
m=meReadMergedImage(mi,1,'si');
ds=mi.imageSize(1)/size(m,1);
[s,freqs]=RadialSpectrum(m,mi.ctf(1),mi.pixA*ds);
mi.ctf(1).B=20;
c2=meGetEffectiveCTF(mi,size(m)).^2;
c1=RadialCTF(c2,mi.ctf(1),mi.pixA*ds);
%
n1=numel(s);
c1dc=mean(s(round(n1*.8):round(n1*.9)));
%%

mxv=.0001;

plot(freqs,[GaussFilt(min(mxv,s),.1) c1*.00006+c1dc]);
title(['151216 ' mi.baseFilename ' B=20'],'interpreter','none');
ylabel('Spectral density of normalized image');
xlabel('Frequency, Å^{-1}');

