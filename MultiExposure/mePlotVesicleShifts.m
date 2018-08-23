function mePlotVesicleShifts(mi,scaleMag)
% plot the vesicle shifts as vectors.
% scaleMag is the magnification of shifts relative to vesicle positions.
% 
if nargin<2
    scaleMag=100;
end;

nim=size(mi.vesicle.shiftX,2);
scl=mi.pixA/10;  % nm per pixel
str=['b-'; 'g-'; 'r-'];
nim=min(nim,size(str,1));
xs=mi.vesicle.x*scl;  % position in nm
ys=mi.vesicle.y*scl;
shx=mi.vesicle.shiftX*scaleMag*scl;  % drift, in nm/scaleMag
shy=mi.vesicle.shiftY*scaleMag*scl;
npts=numel(xs);
clf;
axis([1,mi.imageSize(1),1,mi.imageSize(2)]*scl);
legendStr=cell(0);
hold on;
for j=1:nim
    for i=1:npts
        plot([xs(i)+shx(i,2) xs(i)+shx(i,j)],[ys(i)+shy(i,2) ys(i)+shy(i,j)],str(j,:));
    end;
    legendStr{j}=num2str(j);
end
hold off;
% legend(legendStr);
xlabel('Location, nm');
title(['Vector length x ' num2str(scaleMag) ' Shifts 1-2 Blue,  2-3 Red']);