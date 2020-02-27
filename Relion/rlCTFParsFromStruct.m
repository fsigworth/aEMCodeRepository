function pars=rlCTFParsFromStruct(s,index,B)
% function pars=rlCTFParsFromStruct(s,index,B)
% From the variables in a star file at line index, create a ctfPars struct
% that can be passed to our CTF() function.
% e.g.
% [~,d]=ReadStarFile('micrographs_ctf.star');
% d=d{1};
% P=rlCTFParsFromStruct(d,1);
% C=CTF(n,P);
%  The optional B value is placed into pars.B, otherwise the default B will
%  be used.
defaultB=40;

pars=struct;
pars.lambda=EWavelength(s.rlnVoltage(index));
pars.defocus=(s.rlnDefocusU(index)+s.rlnDefocusV(index))/2e4;
pars.deltadef=(s.rlnDefocusU(index)-s.rlnDefocusV(index))/2e4;
pars.theta=s.rlnDefocusAngle(index)*pi/180;
pars.alpha=s.rlnAmplitudeContrast(index);
pars.Cs=s.rlnSphericalAberration(index);
pars.pixA=s.rlnDetectorPixelSize(index)/s.rlnMagnification(index)*1e4;
if nargin>2
    pars.B=B;
else
    pars.B=defaultB;
end;