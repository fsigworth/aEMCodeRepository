function W = EllipticalVesicleFromModel(n,aMinor,e,model,org)
%  function W = EllipticalVesicleFromModel(n,a,e,model,org)
% Given a density cross-section vector 'model', construct the density of a
% spherical vesicle in an image of size n, with radius a and center org.
% Let nm be numel(model). Then the density of the center of
% the membrane is taken to be model(nm/2+1) when nm is even, or
% model((nm+1)/2) when nm is odd; this is the same convention as fftshift.
% VesicleFromModel(n,a,ones(d,1),org) gives the same result as
% VesicleDensity(n,a,d,org).  n can be a 2-element vector to make a
% rectangular result.
tol = 1e-3;  % minimum density different to invoke calculation of a shell
maxEllipticity=.8;
if numel(n)<2  % scalar means a square image
    n=[n n];
end;
ctr=ceil((n+1)/2);  % correct also for odd n
if nargin<5
    org=ctr;
end;

org=org(:)';  % make it a row vector
t=numel(model)/2;
if sum(abs(e))>maxEllipticity
    e=maxEllipticity*e/sum(abs(e));  % don't allow it to get too large
end;
aBound=aMinor/sqrt(1-sum(abs(e)));
n1=2*ceil(aBound+numel(model)/2)+2;  % minimim square that bounds the vesicle.
if n1<min(n)  % We'll model the vesicle in a smaller square area
    ctr1=ceil((n1+1)/2);
    fracShift=org-round(org);
    org1=fracShift+ctr1;
else
    org1=org;
    n1=n;
end;
if numel(n1)<2
    n1=[n1 n1];
end;
model=[0;model(:);0];  % pad with zeros
nm=numel(model);
maxDensity=max(abs(model));
mctr=ceil((nm+1)/2);

W1=single(zeros(n1));

if ~any(e(:))  % all ellipticities are zero

    r=Radius(n1,org1);
    rm=r-.5;
    rp=r+.5;
    for i=2:nm
        r0=i+aMinor-mctr-.5;
        b=model(i-1)-model(i);
        if r0>0 && abs(b)>tol*maxDensity
            %         W1=W1+b*SphereDensity(n1,r0,org1);
            W1=W1+b*real(rp.*sqrt(r0^2-rp.^2)+r0^2*atan(rp./(sqrt(r0^2-rp.^2)))...
                -rm.*sqrt(r0^2-rm.^2)-r0^2*atan(rm./(sqrt(r0^2-rm.^2))));
        end;
    end;
    
else  % code for elliptical vesicle
    nt=4*n(1);  % number of theta steps
    nume=numel(e);
    e(nume+1:4)=0;
    [r thetas]=Radius(n1,org1);
    rm=r-.5;
    rp=r+.5;
    itheta=round(1+(1+thetas/pi)*nt);
    dtheta=pi/nt;
    aMajor=sqrt(1/(1-(sum(e.^2))))*aMinor;
%     e0=sqrt(1-(aMinor/aMajor)^2)
    
    for i=2:nm
        
        diffDensity=model(i-1)-model(i);
            t=i-mctr-.5;  % distance from the membrane center
        if aMinor+t>0 && abs(diffDensity)>tol*maxDensity
%             Modify the ellipticity to keep uniform thickness
            es=sqrt(1-(aMinor+t)^2/(aMajor+t)^2)...
                *e/sqrt(sum(e.^2));
            a2=(aMinor+t)^2./(1-es(1)*cos(2*dtheta*itheta)-es(2)*sin(2*dtheta*itheta)...
                -es(3)*cos(3*dtheta*itheta)-es(4)*sin(3*dtheta*itheta));
            
            W1=W1+diffDensity*real(rp.*sqrt(a2-rp.^2)+a2.*atan(rp./(sqrt(a2-rp.^2)))...
                -rm.*sqrt(a2-rm.^2)-a2.*atan(rm./(sqrt(a2-rm.^2))))...
                ./sqrt(a2)*(aMinor+t);
        end;
    end;
end;  % elliptical

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
