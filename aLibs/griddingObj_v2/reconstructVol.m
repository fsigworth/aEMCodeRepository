function vol=reconstructVol(vol,imgs,angles,ctfs,snr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructs volume from imgs. Projection directions are defined by angles
% ctfs can be empty, or one frequency doman fftshifted ctf per image
% snr is a single number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the gridder object
 g=gridder('serial');
%Initialize volume
volGrid=g.setVol(vol);
%Insert Images
disp('Back Projection...');
tic;
g.insertImgs(volGrid,imgs,angles,ctfs);
toc;
%Weiner filter at snr
volGrid=volGrid.weinerFilt(snr);
%Get the Weiner filtered projected volume from the grid
vol=g.getVol(volGrid);
