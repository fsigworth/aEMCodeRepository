function m=deFixColumns(m,cols)
% function m=deFixColumns(m,cols)
% Replace values in dead columns given by the array cols.
% For every x value in cols, replace the values of m with the mean of the
% two adjacent columns.  This won't work for adjacent columns or columns at
% the edge!

for i=1:numel(cols)
    col=cols(i);
    m(col,:)=mean(m(col-1:2:col+1,:));
end;
