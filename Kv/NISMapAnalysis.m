% NISMapAnalysis.m

% Based on NISMapSubtraction, but cleaner.

% 1. process the maps

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia')

% mapName1='NIS_I_Na_resampled.mrc';
mapName1='MapsForSubtraction/210705/cryosparc_P1_J55_010_volume_map_sharp(6)_copy.mrc';
mapName2='MapsForSubtraction/210705/cryosparc_P6_J40_010_volume_map_sharp(1).mrc';
pdbName1='MapsForSubtraction/SR_7_1_21_dimer_copy.pdb';

[m1,s]=ReadMRC(mapName1);
m2=ReadMRC(mapName2);

n1=size(m1,1); % it's 144
% Working (cropped) map size is 96
n=96;
us=2; % upsampling for rotate
vs=2; % further upsampling for plotting
nu=us*n;
nv=vs*nu;
dsv=us*vs;

msk1=fuzzymask(n,3,0.45*n,.1*n);

baseRotDegrees=34; % rotation we do with the original map.
basePhi=baseRotDegrees*pi/180;
m1c=Crop(m1,n);
m1cr=rsRotateImage(msk1.*m1c,baseRotDegrees);
m1uc=Downsample(m1cr,n*us);
    m1vcr=Downsample(m1uc,nv);


P=[33.932 2.7759]; % alignment for the J40 Na map.
    m2c=Crop(m2,n);
    m2cr=rsRotateImage(msk1.*m2c,P(1)); %
    fsh=FourierShift(n,[0 0 P(2)]); % shift only in z
    m2crs=real(ifftn(fsh.*fftn(m2cr)));
m2uc=Downsample(m2crs,n*us);
    m2vcr=Downsample(m2uc,nv);

% Handle the pdb file

p1=pdb2mat(pdbName1);
mrot=EulerMatrix(basePhi,0,0);

    ctr0=n1/2;
    p1Coords=[p1.X;p1.Y;p1.Z]/s.pixA-ctr0;
    p1rCoords=mrot*p1Coords;
    p1r.X=p1rCoords(1,:); % Zero-based atom coordinates, in A
    p1r.Y=p1rCoords(2,:);
    p1r.Z=p1rCoords(3,:);

% locate the ions
[u,a,b]=unique(p1.element);
%  elements of u are  {'C'}    {'I'}    {'N'}    {'Na'}    {'O'}    {'S'}
NaPtrs=find(b==4); % i.e. Na
IPtrs=find(b==2);  % i.e. I
ptrsI=[NaPtrs(2:3); IPtrs(1)];
nIons=numel(ptrsI); % pointers to ions in pdb structures p1 or p1r

ligands= cell(3,4);
ligands(1,:)={'OE1' 'GLN' 72 'B'};
% ligands(1,:)={'CD' 'GLN' 72 'B'};
ligands(2,:)={'OH' 'TYR' 144 'B'};
ligands(3,:)={'OG' 'SER' 416 'B'};
nL=size(ligands,1);
ptrsL=zeros(nL,1);
disp('Ligands:');
ligandIon=1;
ligandLabels=cell(nL,1);
for i=1:nL
    disp(ligands(i,:));
    pt=find(strcmp(ligands{i,1},p1.atomName) & strcmp(ligands{i,2},p1.resName) ...
        & ligands{i,3}==p1.resNum & strcmp(ligands{i,4},p1.chainID));
    if numel(pt)>0
        ptrsL(i)=pt(1);
    else
        disp('...not found.');
    end;
    ligandLabels{i}=[ligands{i,2} num2str(ligands{i,3}) ' ' ligands{i,1}];
end;

% Valuable variables

m1vcr, m2vcr, p1r, ptrsI, ptrsL, ligandLabels, dsv, nv






 ctru=nu/2+1;
    m1ucr=ERotate3(m1uc,[phi theta psi],[1 1 1]*ctru,nnInterp);
    % Further upsapling
    m1vcr=Downsample(m1ucr,nv);
    
    if showSi % show synthetic ions
        musr=ERotate3(muSi,[phi theta psi],[1 1 1]*ctru,nnInterp);
        mvsr=Downsample(musr,nv);
    end;
    %     show the projection with ion positions marked
    figure(1);
    imags(sum(m1vcr,3));
    title(num2str(angs(i,:)*180/pi))
    
    % Crearte the rotated coordinates p1r
    % mrot=RotMatrix2(-psiAng);
    mrot=EulerMatrix(phi+basePhi,theta,psi);
    % mrot(3,3)=1;
    ctr0=n1/2;
    p1Coords=[p1.X;p1.Y;p1.Z]/s.pixA-ctr0;
    p1rCoords=mrot*p1Coords;
    p1r.X=p1rCoords(1,:);
    p1r.Y=p1rCoords(2,:);
    p1r.Z=p1rCoords(3,:);
    
    ctrv=nv/2+1;
    hold on;
    plot(p1r.X(ptrsI)*us*vs+ctrv,p1r.Y(ptrsI)*us*vs+ctrv,'bo');
    hold off;
    