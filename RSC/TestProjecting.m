% TestProjecting.m
% Test the projections given the alpha, beta, gamma RSC angles.
% 

degR=pi/180;   % degree radians

load 3KG2mapsub2.9A.mat
n0=size(map,1);
mbnCenter=13;  % center of the membrane in the 3D map


n1=144;
m=zeros(n1,n1,n1);
xm=n1/2-mbnCenter;  % Rotate about the center of the membrane.
xs=n1/2-n0/2;
m(1+xs:n0+xs,1+xs:n0+xs,1+xm:n0+xm)=map;

figure(2);

ShowSections(m);
ks=3;
comp=gridMakePreComp(n1,ks);  % Make the pre-compensation function (a 1D array)
F3=gridMakePaddedFT(m,'grid',comp);  % get the 3D fft in a form for slicing.
    % Once the 3D FT has been made, you can use it repeatedly to make a projection:
%%
    
    figure(1);
clf;
SetGrayscale;



alpha=90;
beta=90;
gamma=45;

for beta=0:2:180;

angs=rsDegToEuler([alpha beta gamma]);

% phi=gamma-90;
% theta=beta;
% psi=alpha+90;
% 
% 
% angs=degR* [phi theta psi];

n=size(m,1);  % m is an n x n x n real volume. n must be a multiple of 4.
ks=3;
P2=gridExtractPlaneE(F3,angs,ks);  % angs is a 3x1 vector of Euler angles (radians)
img=gridRecoverRealImage(P2);     % get the un-padded projection image

% Show the particle's projection viewed from the north pole.
imacs(img);
title(num2str([alpha beta gamma]))
drawnow;
end