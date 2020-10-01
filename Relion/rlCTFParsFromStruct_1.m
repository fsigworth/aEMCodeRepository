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
pars.lambda=300;
pars.defocus=(s.rlnDefocusU(index)+s.rlnDefocusV(index))/2e4;
pars.deltadef=(s.rlnDefocusU(index)-s.rlnDefocusV(index))/2e4;
pars.theta=s.rlnDefocusAngle(index)*pi/180;
pars.alpha=0.1;
pars.Cs=2.7;
pars.pixA=1.068/10000*1e4;
if nargin>2
    pars.B=B;
else
    pars.B=defaultB;
end;