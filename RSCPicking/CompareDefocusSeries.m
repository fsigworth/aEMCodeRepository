% CompareDefocusSeries.m
% Simulate images of AMPARs to see the effect of different defocus series.

% Get mi file names
miPath='/EMWork/Hideki/121210/Simulation/SimDefocusSeries/';
d=dir(miPath);
d(1:2)=[];  % delete the . and .. entries
%

% Get the map
ds=2;
mapName='/Volumes/TetraData/Structures/AMPAR/3KG2map58.mrc';

[origMap s]=ReadMRC(mapName);
mpixA=s.pixA;
nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
nt=ceil(nt1/8)*8;
[map, finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
map=map*pixA;  % approx amplitude correction (V-A scaling)
membraneOffset=-70/pixA;  % downsampled map by 2.

