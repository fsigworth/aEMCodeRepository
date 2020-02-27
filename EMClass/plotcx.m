function plotcx(x,y)

colors=[1 .5 .5 ; .5 .5 1];

plot(x,imag(y),'-','color',colors(2,:));
hold on;
plot(x,real(y),'-','color',colors(1,:));
hold off;