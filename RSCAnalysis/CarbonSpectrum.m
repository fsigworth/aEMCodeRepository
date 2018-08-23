function s=CarbonSpectrum(freqs,pars)
% Standard spectrum of carbon film.  Scale, multiply by effCTF and add to 1 to get
% model spectrum of an image.
if nargin<2
    pars=[1 .031 3.3];
end;
s=pars(1)./((freqs/pars(2)).^pars(3)+1)+1;
