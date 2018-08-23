% TestBGalCompression.m
cd /mnt/md0/data/HDBetagal/

pars=struct;
pars.snrRatio=1000;  % target noise-to-noise ratio
pars.lfCutoff=.1; % radius of frequency domain not fitted
pars.displayOn=1;
pars.falconMasking=1;
CompressZMovies('data','dataZ3', pars);
