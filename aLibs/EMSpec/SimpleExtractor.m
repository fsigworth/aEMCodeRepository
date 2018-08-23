function stack=SimpleExtractor(m,mi,boxSize)
% function stack=SimpleExtractor(m,mi,boxSize)
% Extract particle images from the micrograph m.  boxSize is optional;
% default is a size corresponding to 200 angstroms.

nomBoxA=200;  % default box size in angstroms

ds=mi.imageSize(1)/size(m,1);
if nargin<3
    boxSize=NextNiceNumber(nomBoxA/(mi.pixA*ds));
end;
picks=mi.particle.picks;
np=size(picks,1);
stack=zeros(boxSize,boxSize,np,'single');

for i=1:np
    loc=picks(i,1:2)/ds;
    stack(:,:,i)=ExtractImage(m,round(loc)+1,boxSize);
end;

