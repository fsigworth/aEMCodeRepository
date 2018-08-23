function [sp2 aliImgs]=spAlignTrans(stack,sp,refs)
% function [sp aliImgs]=spAlignTrans(stack,sp,refs)
% Perform translational multi-reference
% alignment on the stack.  If rotations were already given, the translation
% coordinates are rotated appropriately so that the spTransformStack
% function will yield the overall transformation.
% In the returned structure sp2 the fields trans and class
% are updated.

[n ny nim]=size(stack);
nr=size(refs,1);
ds=sp.boxSize/n;  % stack downsampling factor
if nr ~= n % need to change the size of the refs
    refs=spDownsample(refs,n);
end;
xytfr=MRAlign(stack,refs,1);
xytfr
sp2=sp;

% rotate the translation values
for i=1:nim
    trans=xytfr(i,1:2)';
    trans
    rot=sp.rot(i);  % Get the previous rotation angle in degrees
    c=cosd(rot); s=sind(rot);
    newTrans=[c s; -s c]*trans*ds;
    newTrans
    sp2.trans(i,:)=sp2.trans(i,:)+newTrans';
end;

sp2.class=xytfr(:,5);
sp2.cc=xytfr(:,6);

if nargout>1 % if the user requested the aligned stack
%     aliImgs=TransformImages(stack,xytfr);
    [sp0 aliImgs]=spTransformStack(stack,sp);
end;
