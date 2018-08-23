function templates=reMakeTemplates(vols, templateAngles, useParFor)
% function templates=reMakeTemplates(vols, templateAngles, useParFor)
%  Make templates from a map or stack of 3D maps, according to the RSC angles.
%  maps is assumed to be n x n x n x n1 x n2 x ... in size.  The resulting
%  templates are n x n x nangs x n1 x n2... in size, where nangs is
%  size(templateAngles,1).
% 
if nargin<3
    useParFor=false;
end;
%
sz=size(vols);
n=sz(1);
if numel(sz)<4
    sz(4)=1;
end;
nVols=prod(sz(4:end));
nangs=size(templateAngles,1);
ks=3;
templates=zeros([n n nangs sz(4:end)],'single');

comp=gridMakePreComp(n,ks);  % Make the pre-compensation function (a 1D array)
for iVol=1:nVols
    F3=gridMakePaddedFT(vols(:,:,:,iVol),'grid',comp);  % get the 3D fft in a form for slicing.
    if useParFor
        
        parfor i=1:nangs
            angs=rsDegToEuler(templateAngles(i,:));
            P2=gridExtractPlaneE(F3,angs,ks);  % angs is a 3x1 vector of Euler angles (radians)
            templates(:,:,i,iVol)=gridRecoverRealImage(P2);     % get the un-padded projection image
        end;
    else
        for i=1:nangs
            angs=rsDegToEuler(templateAngles(i,:));
            P2=gridExtractPlaneE(F3,angs,ks);  % angs is a 3x1 vector of Euler angles (radians)
            templates(:,:,i,iVol)=gridRecoverRealImage(P2);     % get the un-padded projection image
        end;
    end;
end;
