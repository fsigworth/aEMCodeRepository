function angles=SphereAngles3(nphi,ntheta,nhemi)
% function [angles inds]=SphereAngles(nphi,ntheta,nhemi)
% Create an array 3 x t of angle values which are ~uniform on a hemisphere,
%   At each theta value,the number
% of computed phi values is approximately nphi*sin(theta)
% Angles are given in degrees.  Returned angles are [phi theta 0]
% 
% Examples
% angles=SphereAngles(72,18);  % 5 degree steps, 813 values
% angles=SphereAngles(36,9);  % 10 degree steps, 199 values
% 
% Modified so that the theta values go from 0:dtheta:1-dtheta/2 x 90 x nhemi.
% fs 10 Mar 2012

if nargin<3
    nhemi=1;
end;
nthetaM1=ntheta-1;  % there will be 1 more than this number of theta steps

angles=[];
t=0;
dtheta=180/(nthetaM1+.5);

for i=1:ntheta*nhemi
    % Assigning theta and phi values:
    theta=(i-1)*dtheta;

    NPhiSteps=round(min(nphi,sind(theta)*nphi+1));
    for j=1:NPhiSteps
        % psi values are (j-1)*NPsiSteps*2*pi for j=1..NPsiSteps
        phi=(j-1)/NPhiSteps*360;
        % nangles=nangles+1;
        t=t+1;
        angles(t,:)=[phi theta 0];
    end
end;
