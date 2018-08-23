function m=ColorHisto(yVals,ny)
% Create an image m of size [numel(yVals) x ny] where the bottom yVals
% of each column are set to the x value of each column.
nx=numel(yVals);
m=zeros(nx,ny,'single');
for i=1:nx
    top=round(min(yVals(i),ny));
    m(i,1:top)=i;
end;
