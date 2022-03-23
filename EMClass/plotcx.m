function plotcx(x,y)
if nargin<2
    y=x;
    if any(size(y)==1)
        y=y(:);
    end;
    x=1:size(y,1);
end;
colors=[.9 .4 .4 ; .5 .5 1];

plot(x,imag(y),'-','color',colors(1,:));
hold on;
plot(x,real(y),'-','color',colors(2,:));
hold off;