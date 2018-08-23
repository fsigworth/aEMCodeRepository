function mi=meCreateMicrographInfoStruct10
% function mi=meCreateMicrographInfoStruct10
% Create a blank MicrographInfo structure.
% This structure is intended to be stored in a file with the Matlab save
% function with the full path and name [mi.infoPath mi.baseFilename 'mi.mat']

% -added particle.type fs 12Sep11
% 

mi.version=10;
mi.baseFilename='';
mi.imagePath='';
mi.procPath='';
mi.infoPath='';

mi.imageFilenames={};
mi.imageSize=[1 1]*4096;
mi.pixA=0;

mi.doses=0;  % Doses in e/A^2
mi.camera='YaleF20US';
mi.quality=uint8(zeros(1,8));
mi.notes='';
mi.ctf=struct([]);
mi.spectra={};
mi.mergeMatrix=[];

mi.noiseModelCode={};
mi.noiseModelPars=[];

mi.vesicleModel=[];  % vector giving cross-section densities
mi.vesicle.x=[];
mi.vesicle.y=[];
mi.vesicle.r=[];
mi.vesicle.s=[];

mi.boxSize=0;
mi.particle.x=[];
mi.particle.y=[];
mi.particle.vesicle=[];
mi.particle.type=[];
mi.particle.quality=[];


% Example of typical values:
% 
% mi.version	    10
% mi.baseFilename   '10sep19a00002'
% mi.imagePath      '/Volumes/raid3/liguo-Leginon/10sep19a/00exp-good/ori_MRC/'
% mi.procPath       '/Volumes/raid3/liguo-Leginon/10sep19a/00exp-good/merge/'
% mi.infoPath       '/Volumes/raid3/liguo-Leginon/10sep19a/00exp-good/mi/'

% mi.imageFilenames {'10sep19a00002en.mrc' '10sep19a00002em.mrc' '10sep19a00002ef.mrc'}
% mi.imageSize      [4096 4096]
% mi.pixA           2.53

% mi.doses          [10 20 20]'   % e/A^2 (nim x 1)
% mi.camera         'YaleF20U4000' % Allows choice of correction parameters, prewhitening
% mi.quality	    [0 0 0 0 0 0 0 0]  % eight uint8 values; 0 means unclassified
% mi.notes          'Good ice, but much contamination'
% mi.ctf(i)              % CTF parameters structs i=1..nim 
% mi.spectra{i}.df       % Frequency step of power spectrum for image i.
% mi.spectra{i}.spectrum % 2D spectral density, f=0 in the middle, as used for CTF fitting.
% mi.mergeMatrix(i)		 % merge affine-transformation matrices, 3 x 3 x nim

% mi.noiseModelCode={};  {'shot=p(1);' 'spect=p(2)./f;'}
                         % This is a cell array of strings.  Evaluate these to get
                         % the noise model as a function of frequency f.
% mi.noiseModelPars=[];  % Parameter vector which is passed to the code as p.

% mi.vesicleModel        [0 .8 1.6 1.6 1.6 .8 0]
                         % a vector of scattering values, in volts relative
%                        % to H2O, sampled in pixA steps. Typical numeric
                         % value inside the membrane is 1.6 * pixA.
% mi.vesicle.x(j)		 % x-coords of the center of vesicle j, in units of
%                        % original micrograph pixels [1...imageSize(1)).
% mi.vesicle.y(j)        % vector of vesicle y-coords, in units of
                         % original micrograph pixels [1..imageSize(2)]
% mi.vesicle.r(j)        % vector of vesicle radii, in original pixels.
% mi.vesicle.s(j)		 % vector of model scalings.  These scale factors
%                        % are nominally equal to 2*doses(1)*sigma, where sigma
%                        % is the electron scattering factor .073 rad/VA at 200 keV.

% mi.boxSize             % box size used in picking, e.g. 80
% mi.particle.x(k)		 % x coord, in original image pixels, range 1..imageSize(1)
% mi.particle.y(k)       % y coord, in original image pixels, range 1..imageSize(2)
% mi.particle.vesicle(k) % index (j) of the vesicle containing this particle
% mi.particle.type(k)    % 0 means background, 1 means particle
% mi.particle.quality(k) % quality of particle image
