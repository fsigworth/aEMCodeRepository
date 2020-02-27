function rsShowVesicleNumbers(m,mi)
% display the image m and draw numbers of the indices of vesicles.
color='y';
fontSize=14;

imags(m);
n=size(m);
ds=mi.imageSize(1)/n(1);
hold on;
nv=numel(mi.vesicle.x);
for i=1:nv
    x=double(mi.vesicle.x(i)/ds+1);
    y=double(mi.vesicle.y(i)/ds+1);
    text(x,y,num2str(i),'fontsize',fontSize,'color',color,...
        'horizontalalignment','center','verticalalignment','middle');
end;
hold off;
