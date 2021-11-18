% NISMapLoader.m

% Based on NISMapSubtraction, but cleaner
% Part 1: load the maps and coordinates.
doSave=0j;
% 1. process the maps

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia')

names.map5='211021/cryosparc_P1_J340_005_volume_map_sharp.mrc';
% names.map4='210927/cryosparc_P1_J322_005_volume_map_sharp.mrc';

names.pdb1='MapsForSubtraction/SR_7_1_21_dimer_copy.pdb';
names.pdb2='MapsForSubtraction/model_fit_P8j67-07-12_ek-coot-9.pdb'; % Fig to perrhenate
names.pdb3='MapsForSubtraction/Perrheate_0806.pdb'; % Fig to perrhenate
% names.pdb4='210927/IODO_SR_092721chainA-coot-0.pdb';
names.pdb4='210927/IODO_SR_100121_chainA_real_space_refined_180-coot-1.pdb';
names.pdb5=
disp('Loading maps');
[m1,s]=ReadMRC(names.map1);
m2=ReadMRC(names.map2);
m3=ReadMRC(names.map3);
m4=ReadMRC(names.map4);


n1=size(m1,1); % it's 144
% Working (cropped) map size is 96
n=96;
us=2; % upsampling for rotate
vs=2; % further upsampling for plotting
nu=us*n;
nv=vs*nu;
ctrV=ceil((nv+1)/2); % center of the upsampled map
dsv=us*vs; % final upsampling factor

figure(9);
msk1=fuzzymask(n,3,0.45*n,.1*n);
disp('Rotating and padding map 1');
baseRotDegrees=34; % rotation we do with the original map.
basePhi=baseRotDegrees*pi/180;
m1c=Crop(m1,n);
m1cr=rsRotateImage(msk1.*m1c,baseRotDegrees);
m1uc=Downsample(m1cr,n*us);
    m1v=Downsample(m1uc,nv);
mysubplot(2,2,1);
imags(sum(m1v,3));
title(names.map1,'interpreter','none');
drawnow;

disp('Rotating and padding map 2');
P=[33.932 2.7759]; % alignment for the J40 Na map.
    m2c=Crop(m2,n);
    m2cr=rsRotateImage(msk1.*m2c,P(1)); %
    fsh=FourierShift(n,[0 0 P(2)]); % shift only in z
    m2crs=real(ifftn(fsh.*fftn(m2cr)));
    m2uc=Downsample(m2crs,n*us);
    m2v=Downsample(m2uc,nv);
mysubplot(2,2,2);
imags(sum(m2v,3));
title(names.map2,'interpreter','none');
drawnow;
%%
disp('Rotating and padding map 3');
P=[66.14 2.23];
    m3c=Crop(m3,n);
    m3cr=rsRotateImage(msk1.*m3c,P(1)); %
    fsh=FourierShift(n,[0 0 P(2)]); % shift only in z
    m3crs=real(ifftn(fsh.*fftn(m3cr)));
    m3uc=Downsample(m3crs,n*us);
    m3v=Downsample(m3uc,nv);
mysubplot(2,2,3);
imags(sum(m3v,3));
title(names.map3,'interpreter','none');
drawnow;
basePhi2=P(1)*pi/180;
baseShift2=P(2);

basePhi3=basePhi2;
baseShift3=baseShift2;

disp('Rotating and padding map 4');
P=[34 0];
    m4c=Crop(m4,n);
    m4cr=rsRotateImage(msk1.*m4c,P(1)); %
    fsh=FourierShift(n,[0 0 P(2)]); % shift only in z
    m4crs=real(ifftn(fsh.*fftn(m4cr)));
    m4uc=Downsample(m4crs,n*us);
    m4v=Downsample(m4uc,nv); % rotated, shifted, upsampled.
mysubplot(2,2,4);
imags(sum(m4v,3));
title(names.map4,'interpreter','none');
drawnow;
basePhi4=P(1)*pi/180;
baseShift4=P(2);


%%
% Handle the pdb file
figure(8);
clf;
disp(['Loading the pdb file ' names.pdb1]);
p1=pdb2mat(names.pdb1);
mrot=EulerMatrix(basePhi,0,0);
    ctr0=n1/2;
    p1Coords=[p1.X;p1.Y;p1.Z]/s.pixA-ctr0;
    p1rCoords=mrot*p1Coords;
    p1r=p1;
    p1r.X=p1rCoords(1,:); % Zero-based atom coordinates, in voxels
    p1r.Y=p1rCoords(2,:);
    p1r.Z=p1rCoords(3,:);
    p1r.flipped=false(size(p1r.X));
plot(p1r.X,p1r.Y,'.','markersize',1);


disp(['Loading the pdb file ' names.pdb2]);
p2=pdb2mat(names.pdb2);
mrot=EulerMatrix(basePhi2,0,0);

    ctr0=n1/2;
    p2Coords=[p2.X;p2.Y;p2.Z]/s.pixA-ctr0;
    p2rCoords=mrot*p2Coords;
    p2r=p2;
    p2r.X=p2rCoords(1,:); % Zero-based atom coordinates
    p2r.Y=p2rCoords(2,:);
    p2r.Z=p2rCoords(3,:)+baseShift2;
    p2r.flipped=false(size(p2r.X));

hold on;
plot(p2r.X,p2r.Y,'r.','markersize',1);
hold off;

disp(['Loading the pdb file ' names.pdb3]);
p3=pdb2mat(names.pdb3);
mrot=EulerMatrix(basePhi3,0,0);

    ctr0=n1/2;
    p3Coords=[p3.X;p3.Y;p3.Z]/s.pixA-ctr0;
    p3rCoords=mrot*p3Coords;
    p3r=p3;
    p3r.X=p3rCoords(1,:); % Zero-based atom coordinates
    p3r.Y=p3rCoords(2,:);
    p3r.Z=p3rCoords(3,:)+baseShift3;
    p3r.flipped=false(size(p3r.X));

hold on;
plot(p3r.X,p3r.Y,'g.','markersize',1);
hold off;

disp(['Loading the pdb file 4: ' names.pdb4]);
p4=pdb2mat(names.pdb4);
mrot=EulerMatrix(basePhi4,0,0);
    ctr0=n1/2; % note zero-origin system
    p4Coords=[p4.X;p4.Y;p4.Z]/s.pixA-ctr0; % convert to pixels rel. to center of map
    p4rCoords=mrot*p4Coords; % rotate the coordinates
    p4r=p4;
    p4r.X=p4rCoords(1,:); % Zero-based atom coordinates, in voxels
    p4r.Y=p4rCoords(2,:);
    p4r.Z=p4rCoords(3,:)+baseShift4;
    p4r.flipped=false(size(p4r.X));
hold on;
plot(p4r.X,p4r.Y,'k.','markersize',1);
hold off;
drawnow;

    figure(10);
    clf;
plot(p4r.X,p4r.Y,'b.','markersize',1);
title('pdb4')

drawnow;


%% locate the ions
disp('Locating ions in pdb 1');
[u,a,b]=unique(p1.element);
%  elements of u are  {'C'}    {'I'}    {'N'}    {'Na'}    {'O'}    {'S'}
NaPtrs=find(b==4); % i.e. Na
IPtrs=find(b==2);  % i.e. I
ptrsI1=[NaPtrs(2:3); IPtrs(1)];
nIons=numel(ptrsI1); % pointers to ions in pdb structures p1 or p1r
figure(8)
hold on;
plot(p1r.X(ptrsI1),p1r.Y(ptrsI1),'bo','markersize',10);
hold off;

disp('..pdb 2');
[u,a,b]=unique(p2.element);
NaPtrs=find(b==3);
RePtrs=find(b==5);
HOHPtrs=find(strcmp(p2.chainID,'B'));

NaPtrs=find(strcmp('Na',p2.element))
RePtrs=find(strcmp('Re',p2.element))
HOHPtrs=find(strcmp('HOH',p2.resName))

ptrsI2=[NaPtrs(2); HOHPtrs; RePtrs]
% We ignore the first NaPtr as it's a mirror of Na2
hold on;
plot(p2r.X(ptrsI2),p2r.Y(ptrsI2),'ro','markersize',10);
hold off;

disp('..pdb 3');
% [u,a,b]=unique(p3.element);
% restore these if necessary
    p3r.X=p3rCoords(1,:); % Zero-based atom coordinates, in A
    p3r.Y=p3rCoords(2,:);
    p3r.Z=p3rCoords(3,:)+baseShift3;
    %
NaPtrs=find(strcmp('Na',p3.element));
NaPtrs=NaPtrs(p3r.X(NaPtrs)>0); % get only the right-hand monomer m2
% NaXs=p3r.X(NaPtrs)
RePtrs=find(strcmp('Re',p3.element));
RePtrs=RePtrs(p3r.Y(RePtrs)>0);
% ReXs=p3r.X(RePtrs)
HOHPtrs=find(strcmp('HOH',p3.resName));

% Flip the ones unique to m1 over to m2
fl=HOHPtrs([1 3]);
p3r.flipped(fl)=true;
p3r.X(fl)=abs(p3r.X(fl));
p3r.Y(fl)=abs(p3r.Y(fl));

HOHPtrs=HOHPtrs(p3r.X(HOHPtrs)>0);
% HOHXs=p3r.X(HOHPtrs)
ptrsI3=[NaPtrs HOHPtrs RePtrs];
Xs=p3r.X(ptrsI3)
Ys=p3r.Y(ptrsI3)
hold on;
plot(p3r.X(ptrsI3),p3r.Y(ptrsI3),'ko','markersize',10);
hold off;

disp('..pdb 4');
% [u,a,b]=unique(p4.element);
% restore these if necessary
    p4r.X=p4rCoords(1,:); % Zero-based atom coordinates, in voxels
    p4r.Y=p4rCoords(2,:);
    p4r.Z=p4rCoords(3,:)+baseShift3;
    %
NaPtrs=find(strcmp('Na',p4.element));
IPtrs=find(strcmp('IOD',p4.resName));
% HOHPtrs=find(strcmp('HOH',p4.resName));

ptrsI4=[IPtrs NaPtrs];
nL=numel(ptrsI4);

ligandLabels4=cell(nL,1);
for i=1:nL
% ligands(i,:)=[p4.element(ptrsI4(i))  num2str(p4.atomNum(ptrsI4(i)))];
txt=[p4.element{ptrsI4(i)} ' ' num2str(p4.atomNum(ptrsI4(i)))];
disp(txt);
ligandLabels4{i}=txt;
end;

Xs=p4r.X(ptrsI4)
Ys=p4r.Y(ptrsI4)

figure(10)
hold on;
for i=1:nL
    x=Xs(i);
    y=Ys(i);
    plot(x,y,'ko','markersize',10);
    text(x+.5,y,ligandLabels4{i});
end;
hold off;
drawnow;
figure(9);
hold on;
plot(Xs*dsv+ctrV,Ys*dsv+ctrV,'yo','markersize', 10)
hold off;
drawnow;
% ligands=cell(nL,2);


%% locate the ligands in p1

disp('Finding the ligands in p1...');
ligands= cell(3,4);
ligands(1,:)={'OE1' 'GLN' 72 'B'};
ligands(2,:)={'CD' 'GLN' 72 'B'};
ligands(3,:)={'OH' 'TYR' 144 'B'};
ligands(4,:)={'OG' 'SER' 416 'B'};
ligandIons1=[1 1 1 1];
nL=size(ligands,1);
ptrsL1=zeros(nL,1);
ligandLabels1=cell(nL,1);
for i=1:nL
%     disp(ligands(i,:));
    pt=find(strcmp(ligands{i,1},p1.atomName) & strcmp(ligands{i,2},p1.resName) ...
        & ligands{i,3}==p1.resNum & strcmp(ligands{i,4},p1.chainID));
    if numel(pt)>0
        ptrsL1(i)=pt(1);
    else
        disp([' ' ligands(i,:) '  not found.']);
    end;
    ligandLabels1{i}=[ligands{i,2} num2str(ligands{i,3}) ' ' ligands{i,1}];
    disp(['  ' ligandLabels1{i}]);
end;
figure(8)
hold on;
plot(p1r.X(ptrsL1),p1r.Y(ptrsL1),'kd','markersize',10);
hold off;

%%
% Valuable variables
if doSave
disp('Saving variables...');
save NISMapData.mat m1v m2v m3v m4v p1r p2r p3r p4r ptrsI1 ptrsI2 ptrsI3 ptrsI4 ptrsL1 ligandLabels1 ligandLabels4 ligandIons1 dsv nv s names
disp('done.');
end;
return

