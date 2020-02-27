function vlSet1080Figure(figNo,scale,height)
% function vlSet1080Figure(figNo,scale,height)
if nargin<3
    height=1080;
end;
if nargin<2
    scale=1;
end;
if nargin<1
    figNo=1;
end;

figure(figNo);
h=gcf;

b1=scale*60; % left/lower border
b2=scale*30; % other border
bs=b1+b2;
nx=scale*1920-bs;
ny=scale*height-bs;
pos=h.Position;
h.Position=[pos(1:2) nx+bs ny+bs];
