function imgs=forwardProjectWithAngles(vol,angles,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward projects the vol at pojection angles corresponding to the 
% the Euler angles in angles. This is written to be compatible with
% the old mex_forward_project routine. The only difference is that 
% this uses angles instead of axes
% mode is either 'serial' or 'parallel'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Object Oriented Gridding Starts here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the gridder object
 g=gridder(mode);
%Set the volume into the grid
volGrid=g.setVol(vol);
%Extract images
disp('Forward Projection ...');
tic;
imgs=g.extractImgs(volGrid,angles);
toc;

