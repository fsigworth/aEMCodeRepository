function angs=hrMakeProjectionAngles(angSteps,symmetry)
% function angs=hrMakeProjectionAngles(angSteps,symmetry)
% Make a set of angle triplets in degrees, suitable for rlMakeTemplates()
% given angSteps [phi,theta,psi] in degrees.

if numel(angSteps)==1
    angSteps=angSteps*[1 1 1];
end;

nTheta=round(180/angSteps(2))+1;
dTheta=180/(nTheta-1);
nPsi=ceil(360/angSteps(3));
dPsi=360/nPsi;
psis=(0:dPsi:360-dPsi)';
maxNPhi=ceil(360/(angSteps(1)*symmetry))';

nAngs=0;
angs=zeros(0,3);
for i=1:nTheta
    theta=mod(90+(i-1)*dTheta,180); % start at 90 and go up
    nPhi=ceil(sind(theta)*maxNPhi); % Sample phi sparsely when sin(theta) is small.
    dPhi=360/(nPhi*symmetry);
    for j=1:nPhi
        phi=(j-1)*dPhi;
        % enter all the psi values
        angs(nAngs+1:nAngs+nPsi,:)=[phi*ones(nPsi,1) theta*ones(nPsi,1) psis];
        nAngs=nAngs+nPsi;
    end;
end;
