function E=EulerMatrix(angles)
% function E=EulerMatrix(angles)
% Returns the transformation matrix for rotation of a particle through
% the angles theta, phi and psi.  We follow the definition of 
% B. Heymann et al., JSB 151:196-207, 2005.
% fs 10 Feb 07
% Use:
%	[ newx			[ oldx
%	  newy	= E * 	  oldy
%	  newz ]		  oldz ]
%

phi=angles(1);
theta=angles(2);
psi=angles(3);

% Here are the rotations, in reverse order of their application.
r1=[cos(psi) sin(psi)  0		% Rotate about the z' axis by psi
	-sin(psi)  cos(psi)  0
	0			0			1];

r2=[cos(theta)	0		-sin(theta)	% Rotate about y' by theta
	0			1           0
	sin(theta)  0       cos(theta)];

r3=	[cos(phi)	sin(phi)	0	% About z by phi.
	-sin(phi)	cos(phi)	0
	0			0			1];

E=r1*r2*r3;
