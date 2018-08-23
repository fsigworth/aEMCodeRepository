function [reconVol normVol]=rsDoReconstructionSymmetry(clsSums,clsNorms,angs,symmetry)

[n n1 nrefs]=size(clsSums);
ks=3;
comp=gridMakePreComp(n,ks);

angles=zeros(nrefs,3,symmetry);
 gammaStep=360/symmetry;
degAngles=angs;
for i=1:symmetry
    angles(:,:,i)=rsDegToEuler(degAngles);
    degAngles(:,3)=degAngles(:,3)+gammaStep;
end;

normVol=gridMakeNullFT(n,3);
reconVol=normVol;

% We do Fourier insertions into two volumes.  vol0 just has constant planes
% inserted, and is used for normalization.  vol1 is the actual
% reconstruction volume.
disp('Normalization reconstr');
tic
for i=1:nrefs
    nslice=gridMakePaddedFT(clsNorms(:,:,i));
    for j=1:symmetry
        normVol=gridInsertPlane(nslice,normVol,angles(i,:,j));
    end;
end
toc

disp('Actual reconstr');
for i=1:nrefs
    fslice=gridMakePaddedFT(clsSums(:,:,i));
    for j=1:symmetry
        reconVol=gridInsertPlane(fslice,reconVol,angles(i,:,j));
    end;
end;
