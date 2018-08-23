function D=DistanceMatrix(rs)
% rs is an nx2 matrix of x,y coordinates
% Create the distance matrix D(i,j)=||rs(i,:)-rs(j,:)||

n=size(rs,1);
D=zeros(n,n);
for j=1:n
    ds=sqrt((rs(:,1)-rs(j,1)).^2+(rs(:,2)-rs(j,2)).^2);
    D(:,j)=ds;
end;
