function [projs,angs,shifts]=rlMakeRelionProjections(map,d,indices,shiftPixA,dotCount);
% function [projs,angs,shifts]=rlMakeRelionProjections(map,d,indices,shiftPixA,dotCount);
% Given the struct d from reading a Relion data.star file, compute 2D
% projections from the 3D map, with indices being the line numbers in the
% file. shiftPixA, if nonzero, gives the pixel size of the images for
% shifts determined by the d.rlnOriginAngst variables. Set shiftPixA=0 to
% make no shifts. Dots are written to the command window to show progress.
% dotCount setsthe number of projections made for each dot printed; default
% is 0 which means no dots.
if nargin<3
    indices=[];
end;
if nargin<4
    shiftPixA=0; % i.e. do no shifting
end;
if nargin<5
    dotCount=0; % i.e. don't show dots while running.
end;

% this is what we determined from rlMakeFakedataset
% rl_rot = -phi = -angs(i,1)
% rl_tilt=theta= angs(i,2)
% rl_psi = -psi-90 = -angs(i,3)-90;
% we apply shifts after rotating and projecting.
% But actually we assign the angles one-to-one after all.

np=numel(indices); % use all entries of indices is empty or too large.
if np>numel(d.rlnAngleRot)
    np=numel(d.rlnAngleRot);
elseif np==0
    np=numel(d.rlnAngleRot);
    indices=1:np;
end;
angs=zeros(np,3);
shifts=zeros(np,2);
kShift=0;
if shiftPixA>0
    kShift=1/shiftPixA;
end;

for i=1:np
    j=indices(i);
    angs(i,1)=d.rlnAngleRot(j);
    angs(i,2)=d.rlnAngleTilt(j);
    angs(i,3)=d.rlnAnglePsi(j);
    shifts(i,1)=d.rlnOriginXAngst(j);
    shifts(i,2)=d.rlnOriginYAngst(j);
end;
projs=rlMakeTemplates(angs,map,dotCount,shifts*kShift);
