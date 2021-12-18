function a=LinLeastSquares(F,y,suppressWarnings)
% function a=LinLeastSquares(F,y,suppressWarnings)
% Given a matrix F of function values (each column is values of 1 basis function)
% find the column vector of coefficients a that minimizes
% ||y-Fa||^2
% Thus each column of F corresponds to one of the functions to be fitted.
% The fitted function is obtain as F*a.
% If suppressWarnings=0, Matlab will give warnings about singular matrices.
% Default is suppressWarnings=1.

if nargin<3
    suppressWarnings=1;
end;
[n m]=size(F);
if numel(y)~=n
    error('Number of rows in F and y don''t match');
end;
y=y(:);
if suppressWarnings
    warning('off','MATLAB:singularMatrix');
    a=(F'*F)\(F'*y);
    warning('on','MATLAB:singularMatrix');
else
    a=(F'*F)\(F'*y);
end;
% a'

