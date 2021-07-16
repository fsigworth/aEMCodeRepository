% NISMapLoader.m

% Based on NISMapSubtraction, but cleaner
% Part 1: load the maps and coordinates.

% 1. process the maps

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia')

% mapName1='NIS_I_Na_resampled.mrc';
names.map1='MapsForSubtraction/210705/cryosparc_P1_J55_010_volume_map_sharp(6)_copy.mrc';
names.map2='MapsForSubtraction/210705/cryosparc_P6_J40_010_volume_map_sharp(1).mrc';
names.map3='MapsForSubtraction/210705/cryosparc_P8_J67_010_volume_map_sharp(1).mrc'; % perrhenate

names.pdb1='MapsForSubtraction/SR_7_1_21_dimer_copy.pdb';
names.pdb2='MapsForSubtraction/model_fit_P8j67-07-12_ek-coot-9.pdb'; % Fig to perrhenate

disp('Loading maps');
[m1,s]=ReadMRC(names.map1);
m2=ReadMRC(names.map2);
m3=ReadMRC(names.map3);

n1=size(m1,1); % it's 144
% Working (cropped) map size is 96
n=96;
us=2; % upsampling for rotate
vs=2; % further upsampling for plotting
nu=us*n;
nv=vs*nu;
dsv=us*vs;

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

%%
% Handle the pdb file

disp(['Loading the pdb file' names.pdb1]);
p1=pdb2mat(names.pdb1);
mrot=EulerMatrix(basePhi,0,0);
    ctr0=n1/2;
    p1Coords=[p1.X;p1.Y;p1.Z]/s.pixA-ctr0;
    p1rCoords=mrot*p1Coords;
    p1r=p1;
    p1r.X=p1rCoords(1,:); % Zero-based atom coordinates, in A
    p1r.Y=p1rCoords(2,:);
    p1r.Z=p1rCoords(3,:);
mysubplot(2,2,4);
plot(p1r.X,p1r.Y,'.','markersize',1);


disp(['Loading the pdb file' names.pdb2]);
p2=pdb2mat(names.pdb2);
mrot=EulerMatrix(basePhi2,0,0);

    ctr0=n1/2;
    p2Coords=[p2.X;p2.Y;p2.Z]/s.pixA-ctr0;
    p2rCoords=mrot*p2Coords;
    p2r=p2;
    p2r.X=p2rCoords(1,:); % Zero-based atom coordinates, in A
    p2r.Y=p2rCoords(2,:);
    p2r.Z=p2rCoords(3,:)+baseShift2;
hold on;
plot(p2r.X,p2r.Y,'r.','markersize',1);
hold off;
drawnow;


%% locate the ions
disp('Locating ions in pdb 1');
[u,a,b]=unique(p1.element);
%  elements of u are  {'C'}    {'I'}    {'N'}    {'Na'}    {'O'}    {'S'}
NaPtrs=find(b==4); % i.e. Na
IPtrs=find(b==2);  % i.e. I
ptrsI1=[NaPtrs(2:3); IPtrs(1)];
nIons=numel(ptrsI1); % pointers to ions in pdb structures p1 or p1r
hold on;
plot(p1r.X(ptrsI1),p1r.Y(ptrsI1),'bo','markersize',10);
hold off;

disp('..pdb 2');
[u,a,b]=unique(p2.element);
NaPtrs=find(b==3);
RePtrs=find(b==5);
HOHPtrs=find(strcmp(p2.chainID,'B'));

ptrsI2=[NaPtrs(2); HOHPtrs; RePtrs];
% We ignore the first NaPtr as it's a mirror of Na2
hold on;
plot(p2r.X(ptrsI2),p2r.Y(ptrsI2),'ro','markersize',10);
hold off;



%% locate the ligands
disp('Finding the ligands...');
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
hold on;
plot(p1r.X(ptrsL1),p1r.Y(ptrsL1),'kd','markersize',10);
hold off;



% Valuable variables
disp('Saving variables...');
save NISMapData.mat m1v m2v m3v p1r p2r ptrsI1 ptrsL1 ligandLabels1 ligandIons1 dsv nv s names
disp('done.');
return

