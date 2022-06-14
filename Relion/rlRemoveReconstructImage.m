% rlRemoveReconstructImage.m
% Deletes the rlnReconstructImageName field from the input particles.star
% file, and writes out the result as particles_subOnly.star back into the
% same directory.

inStarName='Extract/job270/particles.star'; % Subtracted stack!

[pa, nm, ex]=fileparts(inStarName);
outStarName=[pa filesep nm '_subOnly' ex];
% Get the pair of micrograph.star files. We assume that the unsub and
% subtracted micrograph names are on the same rows in the files.

if ~exist(inStarName,'file')
    disp(['The input file ' inStarName ' can''t be found.']);
    return
end;

% Get the input particle.star file.
disp(['Reading the input file: ' inStarName]);
[nmp,datp]=ReadStarFile(inStarName);
dOut=datp{2}; % Copy the particle data.

if isfield(dOut,'rlnReconstructImageName')
    dOut=rmfield(dOut,'rlnReconstructImageName');
    disp(['Writing the output file ' outStarName]);
    datp{2}=dOut;
    WriteStarFile(nmp,datp,outStarName);
    disp('done.')
else
    disp('Already ok: the file')
    disp(inStarName)
    disp('has no rlnReconstructImage field. Nothing written.');
end;
