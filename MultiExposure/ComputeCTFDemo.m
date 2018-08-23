% ComputeCTFDemo.m

[name pa]=uigetfile('*','Select a micrograph');

cd(pa);

[m pixA]=ReadEMFile(name);  % This reads TIFF, DM3, MRC, Imagic

if pixA==0  % Didn't get pixel size from the file
    pixA=input('Angstroms per pixel? ');
end;
defocus=input('Approximate defocus, um? ');

m=RemoveOutliers(m);  % if this is a ccd image, remove outlying pixels

% note that the ctffit program assumes square images with even numbers of
% pixels.
ctfPars=meFitCTFs(m,pixA,defocus,0);

ctfPars

