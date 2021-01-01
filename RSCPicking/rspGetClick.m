function [coords, b]=rspGetClick(dis)
% Get a click or keypress and return the coordinates in units of the padded
% micrograph pixels, zero-based.
pointers={'thicksquare' 'squareshadow' };
ptrIndex=dis.readOnlyMode+1;  % squareShadow pointer in normal mode.
% [loc,b]=WaitForClick;
% x=loc(1);
% y=loc(2);
[x, y, b]=Myginput(1,pointers{ptrIndex});
% [x y b]=ginput(1);
%             rawCoords=[x y]

coords=zeros(1,3);
if numel(x)<1 || numel(b)<1
    b=0;
    return
end;
% We'll return the coordinates as absolute micrograph 0-based.
coords(1)=(x+dis.org(1))*dis.ds;
coords(2)=(dis.size(2)-y+dis.org(2)-1)*dis.ds;
% cropOffset=floor((mi.padImageSize-mi.imageSize)/2);
% coords=coords-cropOffset;
% coords
