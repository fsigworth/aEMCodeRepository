% rlMakeCTFInvJpegs.m

starName='CtfFind/job005/micrographs_ctf.star'

pars=struct;
pars.writeMiFile=0; % Write out each mi.txt file.
pars.writeMergedImage=0;
pars.writeMergedSmall=0;
pars.writeJpeg=0;
pars.writeJpegInv=-1; % reverse contrast
pars.compFraction=0.4;
pars.dsSmall=4; % downscaling for Small and jpeg
pars.disFc=0.3;
pars.disHP=.005;
pars.lastLine=inf;
pars.firstPeakAmp=1;

rlStarToMiFiles(starName,pars);

disp('done.')
