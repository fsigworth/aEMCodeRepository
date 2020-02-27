function ShowSections2( map, ctr, yAngle )
% function ShowSections( map, ctr, yAngle);
% Show 3  sections and projections of the 3d map.
% The optional ctr vector gives the [x y z] values at which sections are
% taken.  Default is FFT center, also taken when ctr=[].
% The optional yAngle is the angle in degrees of the y-projection (Y' axis)
% relative to the X-axis; default is 90.
%
SetGrayscale;
[nx ny nz]=size(map);
ctr0=ceil(([nx ny nz]+1)/2);

if nargin<2 || numel(ctr)<3
    ctr=ctr0;
end;

if nargin<3
    yAngle=90;
end;

xs=1:nx;
ys=1:ny;
zs=1:nz;

cx=ctr(1);
cy=ctr(2);
cz=ctr(3);

cx0=ctr0(1);
cy0=ctr0(2);
cz0=ctr0(3);

if yAngle>=90
    yxs=[cx cx];
    yxs0=[cx0 cx0];
    yys=[1 ny];
    yys0=yys;
    rotMap=map;
    yString='Y';
else
    tanY=tand(yAngle);
    yxs=[1 nx];
    yxs0=yxs;
    yys=cy+tanY*[1-cx nx-cx];
    yys0=cy0+tanY*[1-cx0 nx-cx0];
    rotMap=rsRotateImage(map,yAngle-90);
    yString='Y''';
end;

% X-Y plane
subplot(3,3,1);
imacs(xs,ys,map(:,:,cz));
hold on;
plot([1 nx ],[cy cy],'b-');
% plot([cx cx],[1 ny],'g-');
plot(yxs,yys,'g-');

hold off;
xlabel('X');
ylabel('Y');

subplot(3,3,4);
imacs(xs,ys,sum(map,3));
hold on;
plot([1 nx ],[cy0 cy0],'b-');
% plot([cx0 cx0],[1 ny],'g-');
plot(yxs0,yys0,'g-');
hold off;
xlabel('X');
ylabel('Y');

% X-Z plane

subplot(3,3,2);
imacs(xs,zs,map(:,cy,:));
hold on;
plot([1 nx ],[cz cz],'b-');
plot([cx cx],[1 nz],'r-');
hold off;
xlabel('X');
ylabel('Z');
title('Sections');

subplot(3,3,5);
imacs(xs,ys,sum(map,2));
hold on;
plot([1 nx ],[cy0 cy0],'b-');
plot([cx0 cx0],[1 nz],'r-');
hold off;
xlabel('X');
ylabel('Z');
title('Projections');

% Y-Z plane
subplot(3,3,3);
imacs(ys,zs,rotMap(cx,:,:));
hold on;
plot([1 ny ],[cz cz],'g-');
plot([cy cy],[1 nz],'r-');
hold off;
xlabel(yString);
ylabel('Z');

subplot(3,3,6);
imacs(xs,ys,sum(rotMap,1));
hold on;
plot([1 ny ],[cz0 cz0],'g-');
plot([cy0 cy0],[1 nz],'r-');
hold off;
xlabel(yString);
ylabel('Z');

subplot(3,3,7);  % Plot 1D sections along x, y and z
plot(xs,rotMap(:,cy,cz)','b-',...
     ys,rotMap(cx,:,cz),'g-',....
     zs,squeeze(rotMap(cx,cy,:)),'r-');
legend('X','Y','Z');
title('Line sections');

subplot(3,3,8);  % Plot corresponding 1D projections
plot(xs,squeeze(sum(sum(rotMap,3),2)),'b-',...
     ys,squeeze(sum(sum(rotMap,3),1)),'g-',....
     zs,squeeze(sum(sum(rotMap,2),1)),'r-');
title('Line projections');
subplot(3,3,1);  % Go back to here in case the user wants to add a title
