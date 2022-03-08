function plotcx(x,y)
if nargin<2
    y=x;
    x=1:numel(y);
end;
colors=[.9 .4 .4 ; .5 .5 1];

plot(x,imag(y),'-','color',colors(1,:));
hold on;
plot(x,real(y),'-','color',colors(2,:));
hold off;