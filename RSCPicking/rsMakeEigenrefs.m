function eigenSet=rsMakeEigenrefs(templates,nterms)
% function eigenSet=rsMakeEigenrefs(templates,nterms)
% Given a set of templates (nt x nt x nim) where nim is the product of the
% three higher dimensions nHemi, nGamma and nAngs, create an eigenimage
% expansion using nterms.  A typical value for nterms is 20.
% The result is the structure eigenSet with the fields...
% eigenSet.imgs: Eigenimages, nt x nt x nterms, each with unity power.
% eigenSet.vList: Expansion of each template in factor space, nterms x nim.
% eigenSet.vListNorm: the same, but each vector (column) is normalized.
% eigenSet.ampList:  rms amplitude of each vList vector, nim x 1.
%   So that vListNorm = vList / ampList.
% eigenSet.termVar:  Shows the power as a function of the number of terms
%   by computing the cumulative sum of the squared amplitudes, normalized
%   by the original template variances.  nterms x nim
%%
[nt, nt1, nHemi, nGamma, nAngs]=size(templates);  % nHemi should be 2.
nim=nHemi*nGamma*nAngs;

templateVector=double(reshape(templates,nt*nt1,nim));
tplVar=sum(templateVector.^2)';

if nterms>100
    disp('starting full svd');
    [u, s, v]=svd(templateVector,0);  % full svd
    u=u(:,1:nterms);
    s=s(1:nterms,1:nterms);
    v=v(:,1:nterms);
else
    disp('starting partial svd..');
    [u, s, v]=svds(templateVector,nterms);  % takes 5x less time for 80 than for all 4096.
end;
disp('...done.');
% templateVector is npix x nrefs
% u is the npix x nterms matrix of eigenimages (each np pixels image is a column)
% s is the nterms x nterms diagonal matrix of singular values
% v is the nrefs x nterms matrix of coefficients
% so X ~= u.s.v'
%%
% scale the projection vector
vp=v*s';  % Get the absolute scaling for each projection. u*vp gives the templates back.
vpvar=sum(vp'.^2)';  % Get the sum of each row squared (=var of each projection).
projamps=sqrt(vpvar);
vpn=vp./repmat(projamps,1,nterms);  % Normalized projection coordinates.
% the sum of squares of each row of vpn is 1

% Returned quantities
% Eigenimages, nt x nt x nterms, each with unity power.
eigenSet.imgs=single(reshape(single(u),nt,nt1,nterms));
% vList: expansion of each template in factor space, nterms x nim
eigenSet.vList=reshape(vp',nterms,nHemi,nGamma,nAngs);
% vListNorm: expansion of each template, but each vector (column) is normalized.
eigenSet.vListNorm=reshape(vpn',nterms,nHemi,nGamma,nAngs); % nterms x nim
% ampList: rms amplitude of each vList vector.  vListNorm(:,j) = vList(:,j) / ampList(j).
if numel(projamps)>1
    eigenSet.ampList=reshape(projamps,nHemi,nGamma,nAngs);  % nim elements
else
    eigenSet.ampList=projamps*ones(nHemi,nGamma,nAngs);
end;
% Show the power as a function of the number of terms by computing the
% cumulative sum of the squared amplitudes, normalized by the original
% template variances.
eigenSet.termVar=(cumsum(vp'.^2)'./repmat(tplVar,1,nterms))';  % nterms x nim
