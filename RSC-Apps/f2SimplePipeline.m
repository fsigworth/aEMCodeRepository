% f2SimplePipeline.m  Script to run the pipeline

pars.overwrite=1;  % set this if you want to be able to start over from scratch
% otherwise, f2CreateInfoFiles will skip (and not include names) of
% existing files.

names=f2CreateInfoFiles(pars);  % this will put up a file selector, where
% you select the movie_frames directory.

% names=f2FindInfoFiles;  % alternative you can use if mi files already exist.

pars.batchMode=1;  % don't put up file selector for subsequent processing,
                   % use the names variable instead.
pars.useParfor=1;  % Run parallel (with no graphics)

f2DriftTracker(names,pars);

pars.initialDefoci=[3 12];  % Starting defocus guesses
MergeImages(names,pars);


% To run everything interactively, set pars.batchMode=0.  Then
% f2CreateInfoFiles will ask for the movie_files directory, while the other
% programs will ask you to select info files.



