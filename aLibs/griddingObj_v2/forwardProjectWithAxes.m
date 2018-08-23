function imgs=forwardProjectWithAxes(vol,axes,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward projects the vol at pojection angles corresponding to old style
% coordaxes. This is written to be compatible with
% the old mex_forward_project routine.
% mode is either 'serial' or 'parallel'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert axes to angles
angles=axesToEuler(axes);

imgs=forwardProjectWithAngles(vol,angles,mode);


