function [coords, b]=rspGetClick(dis)
% Get a click or keypress and return the coordinates relative to the
% original micrograph, zero-based.

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
coords(1)=(x+dis.org(1))*dis.ds;
coords(2)=(dis.size(2)-y+dis.org(2)-1)*dis.ds;
% coords