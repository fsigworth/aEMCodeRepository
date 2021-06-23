function out=Smooth1D(in,nReps)
% function out=Smooth1D(in,nReps)
% performs a (1 2 1) smoothing on the input column vector, preserving the
% endpoints. If nReps is given, it sets the number of cycles this is
% repeated.

if nargin<2
    nReps=1;
end;

n=numel(in);
for i=1:nReps
    out=in;
    out(2:n-1)=0.5*in(2:n-1)+0.25*(in(1:n-2)+in(3:n));
    in=out;
end;
