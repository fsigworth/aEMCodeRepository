% LoadTomo.m
% Load the jpeg sections of a tomogram and display.
% The display commands are:
%  left button: pick Y
%  u: increase Z
%  d: decrease Z
%  p: show projection
%  q: quit


[fname, pa]=uigetfile('*.jpg','Select tomogram jpegs','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;
%%
if isnumeric(pa)  % user clicked Cancel
    return
end;
cd(pa);

nFiles=numel(fname)/2;
for i=1:nFiles
    m=imread(fname{i});
    [nx,ny,nim]=size(m);
    if i==1
        ms=zeros([nx ny nFiles],'single');
    end;
    ms(:,:,i)=m(:,:,2);
end;
%%
figure(1);
subplot(2,1,1);
SetGrayscale;
imovie(ms);

%%
subplot(1,2,1);
zmin=round(nFiles/4);
zmax=round(nFiles*3/4);
imacs(sum(ms(:,:,zmin:zmax),3));
title('q to exit');
b=0;
z0=round(nFiles/2);
y0=round(ny/2);
zw=1;
while b~='q'
    [x0, y1, b]=ginput(1);
    switch b
        case 'u'
            z0=z0+zw*2;
                mdis=sum(ms(:,:,z0-zw:z0+zw),3);

        case 'd'
            z0=z0-zw*2;
                mdis=sum(ms(:,:,z0-zw:z0+zw),3);
        case 'p'
            mdis=sum(ms(:,:,zmin:zmax),3);
        case 1
            y0=max(1,min(round(y1),ny));
            mdis(:,y0)=max(mdis(:));
    end;
    subplot(1,2,1);
    z0=max(zw+1,min(z0,nFiles-zw));
%     mdis=sum(ms(:,:,z0-zw:z0+zw),3);
    imacs(mdis);
    title(z0);
    zNum=7;
    dz=floor((zNum-1)/2);
    for i=1:zNum
        subplot(zNum,2,2*i);
        imacs(ms(:,y0+i-dz,:));
        axis off;
    end;
end;