function [lambda, sigma, vRel]=EWavelength(kV)
% function [lambda, sigma, vRel]=EWavelength(kV)
% Compute the electron wavelength lambda (in angstroms) 
% the interaction parameter sigma (in radians/V.angstrom) and the 
% relative velocity vRel=v/c,
% given the electron energy kV in kilovolts.
% Uses the relativistic formula 4.3.1.27 in International Tables.  The
% interaction parameter is from eqn. (5.6) of Kirkland 2010.
% The value of c is 299,792,458 m/s.

lambda=12.2639./sqrt(kV*1000+0.97845*kV.^2);

if nargout>1
    u0=510.999;  % electron rest energy in kilovolts
    sigma=2*pi/(lambda*kV*1000)*(u0+kV)/(2*u0+kV);
end;
if nargout>2
    vRel=sqrt(1-(u0./(kV+u0))^2);
end;