function [v,wi]=VesicleFromRings(n,exPos,exSigma,r,org,exAmps,doCrossSection)
% function [v,wi]=VesicleFromRings(n,exPos,exSigma,r,org,exAmps,doCrossSection)
% Compute the 'ring' vesicle model based on the complex radius expansion
% like VesicleFromModelGeneral, but consisting only of the sum of nxc rings
% of positions (in pixels, relative to nominal radius expansion r) if the 1 x nxc array
% exPos, having Gaussian density of std exSigma (1 x nxc).  org is the
% center and exAmps is an nTerms x nxc matrix with cols being the
% amplitudes of the angular expansion of the weights used in summing the
% nxc rings.
% Option: if exAmps is not given, v is returned as a 3D array, n x n x nxc,
% of ring images, and wi is the n^2 x nTerms array of complex exponentials
% This is useful in l.s. fitting to determine exAmps.  The jth weighted
% ring is then given by reshape(v(:,:,j),n2,1).*real(wi*exAmps(:,j)).

if nargin<6 % no amps given, return an n x n x nxc array of ves models.
    returnRingArray=true;
    exAmps=0;
else
    returnRingArray=false;
end;
if numel(n)<2
    n(2)=n;  % assume image is square
end;
v=zeros(n,'single');  % default return
wi=0;
% exAmps=squeeze(exAmps);  % make sure it's a matrix
nTerms=size(exAmps,1);
nxc=numel(exPos);  % number of extra components
if nxc<1
    return
end;
if ~returnRingArray && ~any([1 nxc] == size(exAmps,2))
    warning(['Inconsistent num rows of exAmps vs exPos, ' num2str(size(exAmps)) ' vs ' num2str(numel(exPos))]);
    return
end;

% Compute n1, the size of the working region
rExt=1.2*r;  % increased values for computing n1
rExt(1)=r(1);
org=org(:)';  % make it a row vector
nModel=ceil(max(exPos)-min(exPos)+6*exSigma);  % extent of model including peaks
n1=[1 1]*2*ceil( max(sum(abs(rExt)))+nModel*.7 )+2;  % minimim square that bounds the vesicle.
if n1<min(n)  % We'll model the vesicle in a smaller square area
    ctr1=ceil((n1+1)/2);
    fracShift=org-round(org);
    org1=fracShift+ctr1; 
else
    org1=org;
    n1=n;
end;

n2=prod(n1);  % total pixels to compute

vArray=zeros(n2,nxc,'single');  % 1d array to receive ring models
 nvdx=2*ceil(max(exPos)+4*exSigma)+1; % size of our temp membrane model
% % nvdx=2*ceil(2*exSigma)+1; % size of our temp membrane model
cvdx=ceil((nvdx+1)/2);  % center of membrane model
% Make the 1D profiles and the shell images.
for i=1:nxc
     vdx=Gaussian(nvdx,1,exSigma,cvdx+exPos(i));
% %     vdx=Gaussian(nvdx,1,exSigma,cvdx);
    r1=r;
% %     r1(1)=r1(1)+exPos(i);
    vArray(:,i)=reshape(VesicleFromModelGeneral(n1,r1,vdx,org1,0,doCrossSection),n2,1);
end;

if (returnRingArray && nargout>1) || ~returnRingArray  % compute the complex exponential
        % Initialize the complex exponentials
    [~,theta2D]=RadiusNorm(n1,org1);
    % Work with 1d variables
    theta=reshape(theta2D,n2,1);
    w=exp(1i*theta);  % complex exponential
    wi=ones([n2 nTerms],'single');
    % Get powers of the complex exponential
    for j=2:nTerms
        wi(:,j)=wi(:,j-1).*w(:);
    end;
end;
if returnRingArray
    v1=reshape(vArray,[n1 nxc]);
    v=zeros([n nxc],'single');
    for i=1:nxc
        v(:,:,i)=ExtractImage(v1(:,:,i),round(org),n,1);  % insert the image into the larger one.
    end;
else
    v1=zeros(n2,1);
    for i=1:nxc
        v1=v1+vArray(:,i).*real(wi*exAmps(:,i));
    end;
    v1=reshape(v1,n1);
    v=ExtractImage(v1,round(org),n,1);
end;
