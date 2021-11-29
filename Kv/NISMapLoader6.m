% NISMapLoader6.m
% Based on NISMapSubtraction, but cleaner, and operates on a single map/pdb
% pair.
%
% This is Part 1: load and process the maps and coordinates. Part 2 is
% NISMapAnalysis5.m

doSave=1 % Write out .mat file with data.

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia')

names.map5='211116/cryosparc_P6_J177_005_APO_map_sharp.mrc';
% names.map5='211116/cryosparc_P8_J137_003_REO_map_sharp.mrc'; % The ReO map
names.pdb5='211116/APO_ter.pdb';
names.data5='NISMapData6_APO.mat'; % output file.
% for REO map we should search for 1x NA and 1x RE

disp('Loading the map');
[m5,s]=ReadMRC(names.map5);

n1=size(m5,1); % it's 144
% Working (cropped) map size is 96
n=96;
us=2; % upsampling for rotate
vs=2; % further upsampling for plotting and peak analysis
nu=us*n;
nv=vs*nu;
dsv=us*vs; % final upsampling factor
ctrV=ceil((nv+1)/2); % center of the upsampled map

figure(1); % Show the map, the model, and ion positions circled.

msk1=fuzzymask(n,3,0.45*n,.1*n); % Avoid any artifacts from rotation.

disp('Rotating and padding map 5');
P=[34 0];  % rotation, z-translation
m5c=Crop(m5,n);
m5cr=rsRotateImage(msk1.*m5c,P(1)); % Rotation in XY plane
fsh=FourierShift(n,[0 0 P(2)]); % shift only in Z
m5crs=real(ifftn(fsh.*fftn(m5cr)));
m5uc=Downsample(m5crs,n*us);
m5v=Downsample(m5uc,nv); % rotated, shifted, upsampled.
mysubplot(121);
imags(sum(m5v,3));
title(names.map5,'interpreter','none');
drawnow;

basePhi5=P(1)*pi/180; % Save transformation for operating on map too.
baseShift5=P(2);


%%

disp(['Loading the pdb file 5: ' names.pdb5]);
p5=pdb2mat(names.pdb5);
mrot=EulerMatrix(basePhi5,0,0);
ctr0=n1/2; % shift to put the origin in the center
p5Coords=[p5.X;p5.Y;p5.Z]/s.pixA-ctr0; % convert to pixels rel. to center of map
p5rCoords=mrot*p5Coords; % rotate the coordinates
p5r=p5;
p5r.X=p5rCoords(1,:); % Zero-based atom coordinates, in voxels
p5r.Y=p5rCoords(2,:);
p5r.Z=p5rCoords(3,:)+baseShift5;
p5r.flipped=false(size(p5r.X)); % in some dimer maps, the mirror ion position is used.

mysubplot(122); % plot the pdb atom locations

plot(p5r.X,p5r.Y,'b.','markersize',1);
title(names.pdb5)

%% locate the ions
disp('Locating the ions in pdb 5');

% This time, the sodiums are called 'HOH' :)
ptrsI5=find(strcmp('HOH',p5.resName));
ptrsI5=ptrsI5(3:end); % Pick the last two!

% NaPtrs=find(strcmp('HOH',p5.element));
% IPtrs=find(strcmp('IOD',p5.resName));

% ptrsI5=[IPtrs NaPtrs];
nL=numel(ptrsI5);

ligandLabels5=cell(nL,1);
for i=1:nL
    txt=[p5.element{ptrsI5(i)} ' ' num2str(p5.atomNum(ptrsI5(i)))];
    disp(txt);
    ligandLabels5{i}=txt;
end;
disp('Ion coords:');
Xs=p5r.X(ptrsI5)
Ys=p5r.Y(ptrsI5)

mysubplot(122);
hold on;
for i=1:nL
    x=Xs(i);
    y=Ys(i);
    plot(x,y,'ko','markersize',10);
    text(x+.5,y,ligandLabels5{i});
end;
hold off;
drawnow;

% Mark the ions on the map
mysubplot(121);
hold on;
plot(Xs*dsv+ctrV,Ys*dsv+ctrV,'yo','markersize', 10)
hold off;
drawnow;
% ligands=cell(nL,2);



%%
% Valuable variables
if doSave
    disp(['Saving variables in ' names.data5]);
    save names.data5 m5v p5r ptrsI5 ligandLabels5 dsv nv s names
    disp('done.');
end;
return

