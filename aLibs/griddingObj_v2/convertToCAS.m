function r=convertToCAS(c)
% function r=convertToCAS(c)
% Convert the complex vector or matrix c to a corresponding real one.
% The dimension(s) of c must be even
r=real(c)+imag(c);
