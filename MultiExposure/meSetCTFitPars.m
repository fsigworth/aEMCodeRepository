function CTFitPars=meSetCTFitPars(defoci,pixA,kV,defocusSpread)
% Initializes the CTFitPars array of structs for use by ctfit2.  Defoci is
% a vector of defocus values (in um).  pixA is the number of angstroms per
% pixel.
% We assume a particular form of the B factor: B=B0+defocus*dB.  We also
% hard-wire the minres and max res.
% kV = 200 by default.
%
if nargin<3
    kV=200;
end;
if nargin<4
    defocusSpread=2;
end;
nim=numel(defoci);

% CT fitting parameters
Pa.lambda=EWavelength(kV);
Pa.defocus= 0;
Pa.deltadef=0;
Pa.theta=pi/180*(0:10:80);
Pa.alpha=.07;
Pa.Cs=2;  % 2 mm
Pa.B=0;
Pa.res=pixA;
CTFitPars=Pa;  % Start it off as a struct.
for i=1:nim
    % Based on the nominal defocus, set up the CTF fitting range
    d=defoci(i);
    if d>0
        d=max(d,.1);
    else
        d=min(d,-.1);
    end;
    sp=defocusSpread;
    Pa.defocus=d/sp:d/16:sp*d;
    dfrange = sqrt(d)/2;
%     dfrange=1;
    Pa.deltadef=-dfrange:dfrange/6:dfrange;
    Pa.B=100+100*d;
    
    CTFitPars(i)=Pa;
end;

