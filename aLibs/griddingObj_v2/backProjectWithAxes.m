function vol=backProjectWithAxes(vol,imgs,axes,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Back projects the imgs at pojection corresponding to old style coordaxes
% This is written to be compatible with
% the old mex_back_project routine. 
% mode is either 'serial' or 'parallel'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert axes to angles
angles=axesToEuler(axes);

vol=backProjectWithAngles(vol,imgs,angles,mode);



