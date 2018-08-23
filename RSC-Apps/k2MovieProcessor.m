% k2MovieProcessor.m

% Automatically process movie files.  After starting Matlab, cd to the base
% directory before starting this script.  
% For example, given the path to the first .tif movie file
% 
% 140523/Experiment/movie_files/sq01/xxx.tif
% 
% you would cd to the base directory 140523/Experiment/
% 
% In one of the directories at the sq01 level a .dm4 gain reference file is
% expected to be found.
% 
% In the end the directory structure under Experiment will look like
% 
% Experiment/Info/
% Experiment/Jpeg/
% Experiment/Merged/
% Experiment/Micrograph/
% Experiment/movie_files/
% 
% If the processing is partially done, you can just run this script again
% and it will skip all the files that have already been processed.

k2CreateInfoFiles
gBatchProcessing=1;
k2FindDefocusJump
k2DriftTracker
MergeImages

