function [classIndices, means]=spClassify(stack,eigs,nclasses)
% function [classIndices means]=spClassify(stack,eigs,nclasses)
% Perform k-means classification of the images in the stack, in the basis of
% the eigenimages.  The stack images are assigned to nclasses different
% classes.
% The returned array classIndices gives the class number assigned to each
% image in the stack.  If desired, class means are returned.

[nx, ny, nim]=size(stack);
[nx1, ny1, neigs]=size(eigs);
if nx1 ~= nx  % eigenimages are a different size
    eigs=spDownsample(eigs,[nx ny]);
end;
stack=reshape(stack,nx*ny,nim);
eigs=reshape(eigs,nx*ny,neigs);

fstack=stack'*eigs;  % nim x neigs, images in factor space.

disp('kmeans...');
tic
classIndices=kmeans(fstack,nclasses,'maxiter',300);
disp('done.');
toc

if nargout>1
    % Compute full-sized means
    stack=reshape(stack,nx,ny,nim);
    means=single(zeros(nx,ny,nclasses));
    for i=1:nclasses
        means(:,:,i)=mean(stack(:,:,classIndices==i),3);
    end;
end;
