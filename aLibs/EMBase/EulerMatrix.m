function E=EulerMatrix(phi, theta, psi)
% function E=EulerMatrix(angles)
% function E=EulerMatrix(phi, theta, psi)
%   here angles =[phi theta psi]
% Returns the Euler matrix as defined by Heymann et al., JSB 151:196-207, 2005.
% CCW rotations about z (by phi in the xy plane), then y' (by theta in the zx plane),
% then z'' (by psi in the x''y'' plane). From the point of view of the projection 
% of an object, phi (=rot in Relion) is the spin about the particle's z axis;
% theta is the tilt away from z axis, and psi is the in-plane rotation of the
% projection.
% %
% Returns the transformation matrix for rotation of a coordinate frame through
% the angles phi, theta and psi:
%	
%	[ newx			[ oldx
%	  newy	= E * 	  oldy
%	  newz ]		  oldz ]
%

if nargin<2  % only one argument
    theta=phi(2);
    psi=phi(3);
    phi=phi(1);
end;

r1=[cos(phi) sin(phi)  0		% Rotate about the Z axis by phi
	-sin(phi)  cos(phi)  0
	0			0			1];

r2=	[cos(theta)	0 -sin(theta)	% About Y by theta.
	    0                   1              0
        sin(theta)	0 cos(theta) ];

r3=[cos(psi)    sin(psi)    0
        -sin(psi)   cos(psi)    0
        0               0                 1];
    
E=r3*r2*r1;
