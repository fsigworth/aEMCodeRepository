function E=EulerMatrixInverse(phi, theta, psi)
% function E=EulerMatrixInverse(angles)
% function E=EulerMatrixInverse(phi, theta, psi)
%   here angles =[phi theta psi]
% Returns the inverse of the Euler matrix defined by Heymann et al., 
% JSB 151:196-207, 2005.  The inverse matrix describes the rotation of an object
% or a projection angle rather than the rotation of the coordinate system.
% %
%	[ newx			[ oldx
%	  newy	= E * 	  oldy
%	  newz ]		  oldz ]
%

if nargin<2  % only one argument
    theta=phi(2);
    psi=phi(3);
    phi=phi(1);
end;
E=EulerMatrix(-psi,-theta,-phi);
