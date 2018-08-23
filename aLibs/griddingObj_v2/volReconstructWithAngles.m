function vol=volReconstructWithAngles(vol,imgs,angles,snr,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recontructs Weiner Filtered vol from imgs at projection angles corresponding to the 
% the Euler angles in angles. This is written to be compatible with
% the old mex_back_project routine. The only difference is that 
% this uses angles instead of axes
% mode is either 'serial' or 'parallel'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Object Oriented Gridding Starts here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the gridder object
g=gridder(mode);
volGrid=g.setVol(vol);
%Insert Images
disp('Back Projection...');
tic;
g.insertImgs(volGrid,imgs,angles);
toc;
%Weiner filter at snr
volGrid=volGrid.weinerFilt(snr);
%Get the Weiner filtered volume from the grid
vol=g.getVol(volGrid);


