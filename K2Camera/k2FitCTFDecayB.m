% FitCTFDecayB.m
% Experiment with K2 images of blank beam and thin carbon.
% 
n=3840;
doControl=2;  % 2 means demo image
imageIndex=1;  % 1 low defocus, 2 high defocus

switch doControl
    case 0
        fnum=5;
        %     fnum=1 is 1um; fnum=5 is 2.5 um
        cd('/Users/fred/EMWork/Hideki/140701/control');
        d=dir('Info/');
        cnm=d(fnum+2).name  %  number 3 is 1um; number 7 is 2.5
        mi1=ReadMiFile(['Info/' cnm]);
        imgName=mi1.imageFilenames{1};
        m10=ReadEMFile([mi1.imagePath imgName]);
        
    case 1
        fnum=9;
        cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1');
        d1=dir('Info/');
        dnm=d1(fnum+12).name
        mi1=ReadMiFile(['Info/' dnm]);
        imgName=mi1.imageFilenames{1};
        m10=ReadEMFile([mi1.imagePath imgName]);
        
    case 2  % Read from CTFDemo
        cd('/Users/fred/EMWork/Hideki/140625/CTFDemo')
        dnm='sq10_009_Jun25_16.26.58mi.txt';
        mi1=ReadMiFile(['Info/' dnm]);
        imgName=mi1.imageFilenames{imageIndex};
        m10=ReadEMFile([mi1.imagePath imgName]);
        
end;


%%
Pa=struct;  % CTF fitting parameters
Pa.lambda=EWavelength(200); % lambda in angstroms
Pa.Cs=2;
if imageIndex==2
    Pa.defocus=7:.2:9;
    fRange=[.05 .12]; % frequency range to fit
else
    Pa.defocus=1.8:.2:3;
    fRange=[.06 .2]; % frequency range to fit
end;
Pa.deltadef=-.2:.1:.2;  % astigmatism range
Pa.theta=0:.5:pi;       % astig angle
Pa.alpha=.02;  % This can also be fitted, for use with phase plate.
Pa.B=20:10:60;  % This is fitted if a vector of values is given.

opts=struct;
opts.fExponent=0;  % Assume flat background
opts.k2Mode=1;     % include spectral dip from double-counting
opts.title=imgName;
nu=1024;  % block size (referred to original image sampling)
ds=2;    % downsampling factor

% Select part of an image
n=2048;
c1=201:200+n;
c2=3841-n:3840;
% m1=m10(c,c2);
% m1=m10(c1,c1);
m1=m10;  % select the whole image


opts.blockSize=nu/ds;
fn=RadiusNorm(opts.blockSize)/ds; % normalized frequency
fw=fn/mi1.pixA;  % actual frequency

% Parameters for carbon model
a0=50;
k4=.034;
a4=2080;
expo=3.5;

mtf2=CCDModelMTF(fn,5).^2;
cModel=mtf2.*((k4^expo)./(fw.^expo+k4.^expo)+a0/a4);
% Normalize it to be unity at lower end of range.
fw1=sectr(fw);
% fg=sqrt(prod(fRange));
fg=.05;
p=find(fw1>fg,1);
if numel(p)>0
    cm1=sectr(cModel);
    cModel=cModel/(cm1(p)*exp(-Pa.B(1)*fw1(p)^2));
else
    warning('Couldn''t normalize cModel');
end;
figure(1);
SetGrayscale;

n=size(m1);
m2=Downsample(m1,n/ds);
opts.carbonModel=cModel;
% [P,cc]=FitCTF(m1,Pa,mi1.pixA,4,15,1,opts);
[P, cc]=FitCTFk2(m2,Pa,mi1.pixA,ds,fRange,opts);
P

return


% Fitting spectral data from Li...Agard 2013
% doseRates=[1.023 2.246 4.0542 8.1623 9.7326 12.87 14.9 20.158];
% bVals=[.0048 .0432 .0991 .2223 .2714 .3631 .4226 .5701];
% fit to b as a funciton of doseRates:
% b=-.0274+.0315*doseRate-9.355e-5*doseRate.^2;
% spectrum = N(1 - b*sinc^2(f*2.584)); with N being total counts in image
% Detection fraction:
% detFractions=.866*exp(-doseRates*5.149/400);  % fraction of counts
% Model for spectrum
% f=(0:n/2-1)'/n;
% r=8;  % dose rate
% b=-.0274+.0315*r-9.355e-5*r^2;
% S=(1-b*sinc(f*2.58).^2);
% plot(f,S);




