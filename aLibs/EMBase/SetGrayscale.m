function SetGrayscale()
r = [0:1/255:1]';
h=gcf;
figure(h);
set(h,'color', [1 1 1] );
colormap( [r r r] );
