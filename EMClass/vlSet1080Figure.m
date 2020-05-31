function vlSet1080Figure(figNo,scale,height)
% function vlSet1080Figure(figNo,scale,height)
% function vlSet1080Figure(figNo,scale,size)
% Create a figure for making a 1080p sized figure window.
%  if height is a scalar, the window height is specified, width is 1920.
% The overall size is then scaled by the factor given.
%  if size is 1x2, it it gives height and width.
if nargin<3
    height=1080;
end;
if nargin<2
    scale=1;
end;
if nargin<1
    figNo=1;
end;
    
if numel(height)<2
    height=[1920 height];
end;

figure(figNo);
h=gcf;

b1=scale*60; % left/lower border
b2=scale*30; % other border
bs=b1+b2;
nx=scale*height(1)-bs;
ny=scale*height(2)-bs;
pos=h.Position;
h.Position=[pos(1:2) nx+bs ny+bs];
