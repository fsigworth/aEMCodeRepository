function [max, center, A]=InterpMax3(m)
% function [max, center, A]=InterpMax3(m)
%
% fit the function
% m = max + A.(r-center)^2
% to the matrix of values m(x,y,z).  The returned values
% center refer to index numbers; thus if the optimum occurs
% at m(1,1,1), then 1,1,1 will be returned as center.
% A is the diagonal of the Hessian

[nx, ny, nz]=size(m);
f=zeros(nx,ny,nz,7);
f(:,:,:,1)=ones(nx,ny,nz);
[f(:,:,:,2), f(:,:,:,4), f(:,:,:,6)]=ndgrid(1:nx,1:ny,1:nz);
f(:,:,:,3)=f(:,:,:,2).^2;
f(:,:,:,5)=f(:,:,:,4).^2;
f(:,:,:,7)=f(:,:,:,6).^2;
f=reshape(f,nx*ny*nz,7);
A=f'*f;
y=f'*m(:);
% A=zeros(7,7);
% for i=1:7
% 	for j=1:7
% 		A(i,j)=sum(sum(sum(f(:,:,:,i).*f(:,:,:,j))));
% 	end;
% 	y(i)=sum(sum(sum(f(:,:,:,i).*m)));
% end;
% b=y/A;

b=y'/A;

% now convert the coefficients into the desired form.
A=[b(3) b(5) b(7)];
center=[-b(2)/(2*A(1)) -b(4)/(2*A(2)) -b(6)/(2*A(3))];
% max=b(1)-sum(A.*center.*center);
max=b(1)-A*(center.^2)';