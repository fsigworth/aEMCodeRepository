% TestObj.m
% My version of Fred's test code
%
clear all;
close all;

%Set the path to the griddingObj directory
path(path,'.\griddingObj_v2');
path(path,'.\EMBase');

n=64;   % number of points
na=100;  % number of angles
% random angles
angles=pi*[2*rand(na,1) rand(na,1) 2*rand(na,1)];

% Make a phantom
w=2;  % sharpness of objects
origVol=fuzzymask(n,3,n*.25,w,n*[.4 .4 .5])+2*fuzzymask(n,3,n*.1,w,n*[.5 .8 .5]);

%Display It
figure(1);
ShowSections(origVol);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object Oriented Gridding Starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the gridder object
g=gridder('serial');
%Set the volume into the grid
volGrid=g.setVol(origVol);
%Extract images
disp('Forward Projection ...');
tic;
imgs=g.extractImgs(volGrid,angles);
toc;
%Display images
figure(2);
montage(reshape(imgs,[n n 1 na]),'DisplayRange',[]);
%Clear the original volume grid and set up 
%a grid for volume with all zeros
clear volGrid;
volGrid=g.setVol(zeros(size(origVol)));
%Insert Images
disp('Back Projection...');
tic;
g.insertImgs(volGrid,imgs,angles);
toc;
%Weiner filter at snr
snr=1e-3;
volGrid=volGrid.weinerFilt(snr);
%Get the back projected volume from the grid
vol=g.getVol(volGrid);
%Display it
figure(3);
ShowSections(vol);
drawnow;

