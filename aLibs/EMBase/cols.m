function varargout=cols(x)
% function [x1,x2...]=cols(x)
% Trivial function to pick up the elements of the row vector x, or the
% columns of the matrix x. Use this for example to split coordinates into x
% and y variables.
nc=size(x,2);
varargout=cell(1,nc);
for i=1:nc
    varargout{i}=x(:,i);
end;
