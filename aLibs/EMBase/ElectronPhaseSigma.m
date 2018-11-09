function sigma=ElectronPhaseSigma(kV)
% function sigma=ElectronPhaseSigma(kV)
% electron phase shift per projected potential.
% sigma is returned in radians/(kV*angstrom)
% From Kirkland 1998, p. 65, eqn(5.6)

mc2=511;  % electron rest mass in keV
sigma=2*pi/(EWavelength(kV)*kV)*(mc2+kV)/(2*mc2+kV);
