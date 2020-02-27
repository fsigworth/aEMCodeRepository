function templates=rsMakeTemplates(templateAngles,map)
n=size(map,1);
if mod(n,8)
%     error('Map size must be a multiple of 8');
end;

nangs=size(templateAngles,1);
ks=3;
templates=single(zeros(n,n,nangs));

comp=gridMakePreComp(n,ks);  % Make the pre-compensation function (a 1D array)

F3=gridMakePaddedFT(map,'grid',comp);  % get the 3D fft in a form for slicing.

for i=1:nangs
    angs=rsDegToEuler(templateAngles(i,:));
    P2=gridExtractPlane(F3,angs,ks);  % angs is a 3x1 vector of Euler angles (radians)
    templates(:,:,i)=gridRecoverRealImage(P2);     % get the un-padded projection image
end;
