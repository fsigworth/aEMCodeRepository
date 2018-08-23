function vol=volReconstructWithAxes(vol,imgs,axes,snr,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recontructs Weiner Filtered vol from imgs at projections corresponding to
% old style coordAxes. This is written to be compatible with
% the old mex_back_project routine. 
% mode is either 'serial' or 'parallel'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert axes to angles
angles=axesToEuler(axes);

vol=volReconstructWithAngles(vol,imgs,angles,snr,mode);


