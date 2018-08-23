% VesicleParallaxShifts.m

angle=15;
si=sind(angle);


% Get 2 info files
    [fname, pa]=uigetfile('*mi.*','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
%%
cd(rootPath);
    mis=ReadMiFile([infoPath fname{1}]);
for i=2:numel(fname);
    mis(i)=ReadMiFile([infoPath fname{i}]);
end;

%%
dx=10;
dy=50;
minDist=400;

v1x=mis(1).vesicle.x;
v1y=mis(1).vesicle.y;
v2x=mis(2).vesicle.x;
v2y=mis(2).vesicle.y;
n1x=numel(v1x);
inds2=zeros(n1x,1);
for i=1:n1x
    dists=hypot(v1x(i)-v2x,v1y(i)-v2y);
    [dmin,ind]=min(dists);
    if dmin<minDist
        inds2(i)=ind;
    end;
end;
inds1=(1:n1x)';
inds1(inds2==0)=[];
inds2(inds2==0)=[];
ninds=numel(inds1);

subplot(1,2,1);
plot(v1x(inds1),v1y(inds1),'ro',v2x(inds2),v2y(inds2),'bo');

subplot(1,2,2);
m=meReadMergedImage(mis(1));
xs=0:mis(1).imageSize(1)-1;
ys=0:mis(1).imageSize(2)-1;
imac(xs,ys,imscale(GaussFilt(m,.1),200)+56);
hold on;
plot(v1x(inds1),v1y(inds1),'ro',v2x(inds2),v2y(inds2),'bo');
for i=1:ninds
    z=(v1y(inds1(i))-v2y(inds2(i)))/si*mis(1).pixA/10;
    r=mis(1).vesicle.r(inds1(i))*mis(1).pixA/10;
    str=['  ' num2str(round(z)) ' (' num2str(round(r)) ')'];
    text(double(v1x(inds1(i))),double(v1y(inds1(i))),str);
end;
hold off;
title('Height (radius) in nm');
xlabel([num2str(mis(1).pixA/10) ' nm per pixel']);

m2=meReadMergedImage(mis(2));
subplot(1,2,1);
imacs(GaussFilt(m,.1)-GaussFilt(m2,.06));
drawnow;
%
m1s=imscale(GaussFilt(m,.1));
m2s=imscale(GaussFilt(m2,.06),220)+30;
dt=.1;
for i=1:100
    imac(m1s);
    drawnow;
    pause(dt);
    imac(m2s);
    drawnow;
    pause(dt);
end;
%%
imacs(GaussFilt(m,.1)-2*GaussFilt(m2,.06));
title([mis(1).baseFilename ',  ' mis(2).baseFilename],'interpreter','none');
