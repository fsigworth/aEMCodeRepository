function [v,wi]=VesicleFromRings(n,exPos,exSigma,r,org,exAmps)
% function [v,wi]=VesicleFromRings(n,exPos,exSigma,r,org,exAmps)
% Compute the 'ring' vesicle model based on the complex radius expansion
% like VesicleFromModelGeneral, but consisting only of nxc rings of positions
% (in pixels, relative to nominal radius) exPos, having Gaussian density of
% std exSigma.  org is the origin and exAmps is a matrix with rows being
% amplitudes of the angular expansion of the weights used in summing the
% rings.
% If exAmps is not given, v is returned as an n x nxc array of rings, and
% wi is the n2 x nTerms array of complex exponentials.  This is useful in 
% l.s. fitting to determine exAmps.  A weighted ring is
% then given by reshape(v(:,j),n2,1).*real(wi*exAmps(i,:)').

if nargin<6 % no amps given, return an n2 x nxc array of ves models.
    returnRingArray=true;
    exAmps=0;
else
    returnRingArray=false;
end;
if numel(n)<2
    n(2)=n;  % assume image is square
end;
n2=prod(n);
exAmps=squeeze(exAmps);  % make sure it's a matrix
nTerms=size(exAmps,1);
nxc=numel(exPos);  % number of extra components
if ~returnRingArray && nxc ~= size(exAmps,2)
    error(['Inconsistent num rows of exAmps vs exPos, ' num2str(size(exAmps)) ' vs ' num2str(numel(exPos))]);
end;

vArray=zeros(n2,nxc,'single');
nvdx=2*ceil(max(exPos)+4*exSigma)+1; % size of our temp membrane model
cvdx=ceil((nvdx+1)/2);  % center of membrane model
for i=1:nxc
    vdx=Gaussian(nvdx,1,exSigma,cvdx+exPos(i));
    vArray(:,i)=reshape(VesicleFromModelGeneral(n,r,vdx,org),n2,1);
end;

if (returnRingArray && nargout>1) || ~returnRingArray  % compute the complex exponential
        % Initialize the complex exponentials
    [~,theta2D]=RadiusNorm(n,org);
    % Work with 1d variables
    theta=reshape(theta2D,prod(n),1);
    w=exp(1i*theta);  % complex exponential
    wi=ones([n2 nTerms],'single');
    % Get powers of the complex exponential
    for j=2:nTerms
        wi(:,j)=wi(:,j-1).*w(:);
    end;
end;
if returnRingArray
    v=reshape(vArray,[n nxc]);
else
    v=zeros(n2,1);
    for i=1:nxc
        v=v+vArray(:,i).*real(wi*exAmps(:,i));
    end;
    v=reshape(v,n);
end;
