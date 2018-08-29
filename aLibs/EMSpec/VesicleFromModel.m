function W = VesicleFromModel(n,a,model,org,term2,doCrossSection)
%  function W = VesicleFromModel(n,a,model,org,term2)
% Given a density cross-section vector 'model', construct the density of a
% spherical vesicle in an image of size n, with radius a and center org.
% Let nm be numel(model). Then the density of the center of
% the membrane is taken to be model(nm/2+1) when nm is even, or
% model((nm+1)/2) when nm is odd; this is the same convention as fftshift.
% VesicleFromModel(n,a,ones(d,1),org) gives the same result as
% VesicleDensity(n,a,d,org).  n can be a 2-element vector to make a
% rectangular result.
% 
% term2 (default 0) is extra curvature near the edge of the vesicle.
% Subtracting models with and without term2 gives a constant amplitude when
% term2 is proportional to 1/sqrt(a).  The curvature is changed by letting 
% r -> r+r^2*term2/a, then changing a to a(1+term2).
% if term2=.1, the curvature at the edge is about 10% greater.
tol = 1e-3;  % minimum density different to invoke calculation of a shell

if numel(n)<2  % scalar means a square image
    n=[n n];
end;
ctr=ceil((n+1)/2);  % correct also for odd n
if nargin<4
    org=ctr;
end;
if nargin<5
    term2=0;
end;
if nargin<6
    doCrossSection=false;
end;

org=org(:)';  % make it a row vector
n1=[1 1]*2*ceil(a+numel(model)/2)+2;  % minimim square that bounds the vesicle.
if n1<min(n)  % We'll model the vesicle in a smaller square area
    ctr1=ceil((n1+1)/2);
    fracShift=org-round(org);
    org1=fracShift+ctr1; 
else
    org1=org;
    n1=n;
end;

model=[0;model(:);0];  % pad with zeros
nm=numel(model);
maxDensity=max(abs(model));
mctr=ceil((nm+1)/2);
W1=single(zeros(n1));

% Initialization code from SphereDensity
r=Radius(n1,org1);

% introduce the nonlinearity r := r+(t2/a)*r^2
if term2~=0
r=r.*(1+term2/a*r);  % increasing curvature at the edges
a=a*(1+term2);  % modify the nominal radius to match the final radius.
dr=(1+0.5*term2)*(1+term2); % slope of r near r=a.
else
    dr=1;
end;
% dr=1;%%%%%

if doCrossSection
    for i=2:nm
    b=model(i-1)-model(i);
    r0=a+(i-mctr-.5)*dr;
    if r0>0 && abs(b)>tol*maxDensity
        W1=W1+b*max(0,min(1,(r0-r)*dr));
    end;
    end;
else  % 3D density projection
    
rm=r-.5*dr;
rp=r+.5*dr;
for i=2:nm
    r0=a+(i-mctr-.5)*dr;
    b=model(i-1)-model(i);
    if r0>0 && abs(b)>tol*maxDensity
%         W1=W1+b*SphereDensity(n1,r0,org1);
      W1=W1+b*real(rp.*sqrt(r0^2-rp.^2)+r0^2*atan(rp./(sqrt(r0^2-rp.^2)))...
                  -rm.*sqrt(r0^2-rm.^2)-r0^2*atan(rm./(sqrt(r0^2-rm.^2))));
    end;
end;

end;


if n1<n
    W=ExtractImage(W1,round(org),n,1);  % insert the image into the larger one.
else
    W=W1;
end;
    
% % This code should be faster, but has some error.
% 
% model=[0;0; model(:)];  % pad with zeros
% nm=numel(model);
% mctr=ceil((nm+1)/2);
% W=zeros(n,n);
% rp=Radius(n,org);
% for i=3:nm
%     r0=i+a-mctr-2;
%     b=model(i-2)-2*model(i-1)+model(i)
%     if r0>0 && b~=0
% W=W-b*real(rp.*sqrt(r0^2-rp.^2)+r0^2*atan(rp./(sqrt(r0^2-rp.^2))));
%     end;
% end;
