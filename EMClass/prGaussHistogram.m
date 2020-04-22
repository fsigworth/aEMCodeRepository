% GaussHistogram.m

sigma=1;
mx=3.5;
w=.2;
ctr=0.6;
mHalf=w*ceil(mx/w);
xs=-mHalf:w:mHalf;
ys=w*(2*pi*sigma^2)^-1/2 * exp(-(xs-ctr).^2/(2*sigma^2));
bar(xs,ys,'edgecolor','w','barwidth',.9);
axis([xs(1) xs(end) 0 max(ys)*1.1]);
% ylabel('P(X_j|R_j)');
% xlabel('X_j');