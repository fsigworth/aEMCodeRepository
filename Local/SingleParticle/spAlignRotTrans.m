function [sp, aliImgs]=spAlignRotTrans(stack,sp,refs)
% function [sp, aliImgs]=spAlignRotTrans(stack,sp,refs)
% Perform simultaneous rotational and translational multi-reference
% alignment on the stack.  The input sp structure is used only to give the
% pixel size.  In the returned structure sp the fields trans, rot, flip,
% class and cc are assigned.

n=size(stack,1);
nr=size(refs,1);
ds=sp.boxSize/n;  % stack downsampling factor
if nr ~= n % need to change the size of the refs
    refs=spDownsample(refs,n);
end;
xytfr=MRAlign(stack,refs,3);
sp.trans=xytfr(:,1:2)*ds;  % translations are in pixels of original stack
sp.rot=xytfr(:,3);  % in degrees
sp.flip=xytfr(:,4);
sp.class=xytfr(:,5);
sp.cc=xytfr(:,6);

if nargout>1 % if the user requested the aligned stack
    aliImgs=spTransformStack(stack,sp);
end;
