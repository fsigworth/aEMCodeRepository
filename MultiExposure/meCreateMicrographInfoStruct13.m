function mi=meCreateMicrographInfoStruct13
% function mi=meCreateMicrographInfoStruct13
% Create a blank MicrographInfo structure.
% This structure is intended to be stored in a file with the Matlab save
% function with the full path and name [mi.basePath mi.infoPath mi.baseFilename 'mi.mat']

% derived from version 12, fs July 2013
% 
% Notes on images:
% 
% Image coordinates start with (0,0) in the lower left corner, with the
% first coordinate being the horizontal position.  This follows the EM
% convention.  Particle and vesicle coordinates are given with this
% convention, with the physical step size in angstroms being mi.pixA
% Normally you can find the downsampling factor ds of a derived image m, you can
% compute ds=mi.imageSize(1)/size(m,1).
%
% jpeg images are stored by WriteJpeg() such that they are displayed in the same
% orientation as a binary image is displayed by imacs().
%
% The merged image has units of e/A^2.  that is, meReadImage scales the raw
% image as m=(rawImage-mean(rawImage))/(cpe*pixA^2) where cpe is counts per
% electron.
% The relative contrast of a merged image is gotten as
%   mNormed=m/mi.doses(1)*shotFactor.
% 
% To get the factor by which the merged images were scaled down in v. 11:
%   sp=ifftshift(CCDModelSpectrum2D(iCamera)); shotFactor=sqrt(sp(1,1));
% 
% To scale the image so that the shot noise density is unity, do this
%   mNoiseNormed=m/sqrt(mi.doses);
% 
% To compare the image with the CTF theory, you can use
% the vesicle fitting.  mi.vesicle.s is the factor by which you multiply
% the vesicle model (in V-angstroms) to match m.  Thus the ratio of
% observed to theoretical contrast is
% mi.vesicle.s/(mi.doses(1)*2*sigma/shotFactor) where sigma = .73e-3 
% radians/V.A at 200 keV.  mi.vesicle.s/mi.doses(1) should be about 3e-3.
% 
% Given a projection P of a model in V.A, the corresponding mNormed will be
% modeled by P x effectiveCTF x 2 x sigma.


mi.version=13;
mi.identifier=rand;
mi.baseFilename='';
mi.basePath='';     % global path
mi.imagePath='';    % e.g. 'Images/'
mi.procPath='';
mi.infoPath='';
mi.tempPath='';
mi.stackPath='';

mi.imageFilenames={};
mi.imageSize=[1 1]*0;  % Original micrograph dimension (x,y)
mi.pixA=0;          % Original micrograph pixel size

mi.doses=0;  % Doses in e/A^2
mi.weights=0;
mi.kV=0;
mi.camera=1;
mi.cpe=0;
mi.quality=uint8(zeros(1,8));
mi.notes='';

mi.ctf=struct([]);
mi.mergeMatrix=[];

mi.mask=struct('merge',[],'encoding',[],'data',[]);

mi.noiseModelCode={};
mi.noiseModelPars=[];

mi.vesicleModel=[];  % vector giving cross-section densities
% mi.vesicleModelSampInt=0;  % angstroms per sample, no longer used.

mi.vesicle.x=[];       % nves x 1 vector
mi.vesicle.y=[];
mi.vesicle.r=[];
mi.vesicle.s=[];
mi.vesicle.ok=[];      % nves x 4 matrix of booleans
%  1: found; 2: in bounds in Vesicle_finding_GUI
%  3: successfully fit by rsRefineVesicleFits2; 4: not deselected in picker.
mi.vesicle.shiftX=[];  % nves x nim. increment to x(i) for image j.
mi.vesicle.shiftY=[];  %

mi.boxSize=0;
mi.particle.picks=[];
mi.particle.autopickPars=[];

mi.log=cell(0); % Cell array of strings: processing log


% Example of typical values:
% 
% mi.version	    12
% mi.identifier     .81472   % unique number
% mi.baseFilename   '10sep19a00002'
% mi.basePath       /Volumes/raid3/liguo-Leginon/10sep19a/00exp-good/
% mi.imagePath      'Micrograph/'
% mi.procPath       'Merge/'
% mi.infoPath       'Info/'
% mi.tempPath       'Temp/'
% mi.stackPath      'Stack/'

% mi.imageFilenames {'10sep19a00002en.mrc' '10sep19a00002em.mrc' '10sep19a00002ef.mrc'}
% mi.imageSize      [4096 4096]
% mi.pixA           2.53
% mi.weights        [1 1 1]      % weighting of merged images
% mi.doses          [10 20 20]   % e/A^2 (nim x 1)
% mi.keV            200
% mi.camera         'YaleF20U4000' % Allows choice of correction parameters, prewhitening
% mi.quality	    [0 0 0 0 0 0 0 0]  % eight uint8 values; 0 means unclassified
% mi.notes          'Good ice, but much contamination'
% mi.beamRadius=0;       % If part of the micrograph is black from misaligned
% mi.beamCenter=[0 0];   % beam, give the radius and position, in pixels.
% mi.ctf(i)              % CTF parameters structs i=1..nim.  For example:
%                          ctf(1).lambda=.025;
%                          ctf(1).defocus=2;
%                          ctf(1).Cs=2;
%                          ctf(1).B=200;
%                          ctf(1).alpha=.07;
%                          ctf(1).deltadef=.1;
%                          ctf(1).theta=1.7;
%                          ctf(1).ampFactor=1;   % scaling factor for
%                             vesicle fits.
% mi.mergeMatrix(i)		 % merge affine-transformation matrices, 3 x 3 x
%                          nim, scaled in resolution-independent fashion.
%                        Example values:
%                             mergeMatrix(:,:,1) =
%                                  1     0     0
%                                  0     1     0
%                                  0     0     1
%                             mergeMatrix(:,:,2) =
%                                    1.0124   -0.0016779   3.4273e-05
%                                 0.0015089       1.0126   -0.0015808
%                                         0            0            1
%                             mergeMatrix(:,:,3) =
%                                    1.0618   -0.0076012   0.00018839
%                                 0.0073995       1.0628    -0.011266
%                                         0            0            1
% mi.noiseModelCode={};  {'shot=p(1);' 'spect=p(2)./f;'}
                         % This is a cell array of strings.  Evaluate these to get
                         % the noise model as a function of frequency f.
% mi.noiseModelPars=[];  % Parameter vector which is passed to the code as p.

% mi.vesicleModel        [0 .8 1.6 1.6 1.6 .8 0]
                         % a vector of scattering values, odd number of elelments,
%                        in volts relative
%                        to H2O, sampled in steps of mi.pixA, that is original
%                        pixel size.
%                        Typical numerical value inside the membrane is 1.6.

% mi.vesicle.x(j)		 % x-coord of the center of vesicle j, in units of
%                        % original micrograph pixels [0...imageSize(1)-1).
% mi.vesicle.y(j)        % vector of vesicle y-coords, in units of
                         % original micrograph pixels [0..imageSize(2)-1]
% mi.vesicle.r(j)        % vector of vesicle radii, in original pixels
% mi.vesicle.s(j)		 % vector of model scalings.  These scale factors
%                        % are nominally equal to 2*doses(1)*sigma, where sigma
%                        % is the electron scattering factor .073 rad/VA at 200 keV.
% mi.vesicle.shiftX(j,i) % increment to x(j) for image i.
% mi.vesicle.shiftY(j,i) %

% mi.boxSize             % box size used in picking, e.g. 80
% mi.particle.picks      % each row: x, y, flag, vesicle index, 0 0 0 0.
%                        % x and y in units of pixA, 0,0 = lower left.
%         1700         352          32         190           0           0           0           0
%         1863         300          32         156           0           0           0           0
% mi.particle.autopickPars  % Parameters used by the SimpleRSPicker.
%                           [0.6 1 3.2 35 150]
% mi.log     {'meRefineVesicleFits 2013-04-06 17:18:23'}
