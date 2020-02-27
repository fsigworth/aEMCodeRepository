function ml=Laplacian(m);
% n-dimensional discrete laplacian

nd=ndims(m);
sz=size(m);

ml=-2*nd*m;
for i=1:nd;
    ml=ml+circshift(m,1,i)...
         +circshift(m,-1,i);
end;
