function [compositeMap, proteinMap]=SolventAndProteinDensity(atomCoords, atomTypes, mode, minBox)
%   [compositeMap, proteinMap=SolventAndProteinDensity;  % put up a file
%      selector for a pdb file.
%   [compositeMap proteinMap]=SolventAndProteinDensity(pdbFilename)
%   [compositeMap proteinMap]=SolventAndProteinDensity(atomCoords, atomTypes)
% [compositeMap, proteinMap]=SolventAndProteinDensity(atomCoords, atomTypes, mode, minBox)
% From coordinates and atom types obtained from a pdb file, create a 3D
% density map compositeMap which including simulation of solvent.  The
% protein-alone density is returned as proteinMap.  We use the ReadPDBAtoms
% to get the coordinates and types if a filename is given.  If atomCoords
% is given but not atomTypes, all atoms are assigned the scattering
% potential of oxygen. The resulting maps have a voxel size of 1A.
% The map values are electrostatic potential in volt-angstroms.
% Use the DownsampleGeneral function to change the voxel size post hoc:
%   compositeMap=compositeMap - compositeMap(1);  % remove constant at boundary.
%   m=DownsampleGeneral(compositeMap,n0,1/pixA)*pixA;
% This converts to a map m that is n0^3 voxels of size pixA. The final
% multiplication by pixA means that a projection (sum along projection
% direction) of this map will have units of volt-angstroms.
% Zhiguo Shang and F. Sigworth 2011
originAtEdge=1; %%%%%%%%%%%
if nargin<4
    minBox=0; % create a centered map.
end;
if nargin<3
    mode=1; % smooth grid function
end;
useStepFunction=~mode;
if nargin<1
    [atomCoords,path]=uigetfile('*.pdb');
    if isnumeric(path)
        return
    end;
    cd(path);
end;
if ischar(atomCoords)  % a single string; must be a filename
    pdbName=atomCoords;
    if FileExists(pdbName)
        [atomCoords, atomTypes]=ReadPDBAtoms(pdbName);
    else
        error(['Can''t read file: ' pdbName]);
    end;
end;
    [nd, na0]=size(atomCoords);

if nd ~=3
    error('atomCoords must be a 3 x n array');
end;
if nargin<2
    atomTypes=repmat('O',na0,1);
end;

gridfactor=1;  % Assume 1A grid spacing. (At present the code can't handle changes in this value.)
% Create the maps for atom radius and scattering strength
elements=                     {'H', 'N','C','O','S','P' };
% Van der Walls radii:
radius=containers.Map(elements,[ 0 1.55 1.7 2.0 1.8 1.8 ]);
% Scattering amplitudes in V.A^3, including mean numbers of
% hydrogens associated with each heavy atom type.
weight=containers.Map(elements,([25  108 130  97 100 267]...
    +[0  1.1 1.3 0.2 0.6  0 ]*25));
% Compute the mean scattering density of water
MolesPerM3=1e6/18*0.97;  % assume ice density is 97% that of water
WatersPerA3=1e-30*MolesPerM3*6.022e23;
% Scattering amplitude of the water continuum per 1A^3 grid cube
WaterDensity=(weight('O')+2*weight('H'))*WatersPerA3;

% Scan the atom types.  We generally use only the first character of the types.
% However we specifically ignore sodium (SOD) and chloride (CLO), and we
% skip all hydrogens.
na1=size(atomTypes,1);  % number of atoms to scan
if na1 < na0
    error('Number of atomTypes elements is too small.');
end;
type=char(na1,1);
j=0;
for i=1:na1
    t=atomTypes(i,1);
    type(i)=t;
    switch t
        case 'C'
            if atomTypes(i,2)=='L'  % a chloride
                type(i)='X';
            end;
        case 'S'
            if atomTypes(i,2)=='O'  % sodium
                type(i)='X';
            end;
        case 'H' % skip all hydrogens too.
            type(i)='X';
        otherwise
            if ~isKey(radius,t)
                type(i)='X';
            end;
    end;
end;
q=(type~='X');           % find all valid entries
a0=atomCoords(:,q);       % coords of valid entries
type=type(q);            % atom types of valid entries
na=numel(type);          % number of atoms

% put protein into a solvent box, which will be created in the next step;
mincoord=min(a0,[],2);  % this is vector of the minimum coordinates
maxcoord=max(a0,[],2);

if minBox % Make the smallest box that will fit, regardless of centering
ctrcoord=(mincoord+maxcoord)/2;
span=max(max(maxcoord-mincoord));  % largest span
else
    ctrcoord=[0 0 0]';
    span=2*max(max(abs(maxcoord),abs(mincoord)));
end;
boundary=15;  % number of A of boundary
n=2*ceil((span+2*boundary)*gridfactor/2);
a=(a0+n/2+1-repmat(ctrcoord,1,na))*gridfactor;  % center the coordinates in our box.

if originAtEdge  % Force the zero coordinate to be at the low corner of the box:
    disp('Origin at edge');
    
span=max(max(maxcoord));  % largest coord
n=2*ceil(span/2+boundary);
a=a0;
end;


%% Obtain the minimium distance map of solvent grid mask and EM map caused
% by protein
asigma=1;  % standard deviation of the atomic density Gaussians

nw2=ceil(7*gridfactor); % Kernel size for computing the minimium distance map;
nw=6; % Kernel size for computing the electron density of each atoms;
vols=10*ones(n,n,n);  % Real-space padded volume
volr=vols; % Radius of nearest protein atom. The distance from this atom to grid i
% represents the minimium distance;
vola=char(zeros(n,n,n)); % Type of closest protein atom.

% Generate nw2 and nw size cores for the calculation of minimium distance and EM map;
cint=floor(a);  % Integer coordinates
cfrac=(a-cint)+nw2;  % coordinate remainders, offset to center of kernel box
cfrac1=(a-cint)+nw;

proteinMap=zeros(n,n,n);% density due to molecule's internal atoms
ka=1/(2*asigma^2);   % multiplier in Gauss exponent
kn=(sqrt(2*pi)*asigma)^-3; % Gauss normalization factor

disp(['Creating protein map, ' num2str(na) ' atoms total']);
for ia=1:na  % loop over the atoms
    if mod(ia,5000)==0
        disp(ia);  % show progress
    end;
    i=cint(1,ia);  % Get integer coordinates
    j=cint(2,ia);
    k=cint(3,ia);
    
    i1=i-nw2+1;    % Get box bounds
    i2=i+nw2;
    j1=j-nw2+1;
    j2=j+nw2;
    k1=k-nw2+1;
    k2=k+nw2;
    
    t=type(ia);
    % Calculate the distance map
    r=Radius3(2*nw2,cfrac(:,ia)); % minimium distance from grid map to atom ia.
    indexm=(r-radius(t)<vols(i1:i2,j1:j2,k1:k2)); % All points closer to the present atom.
    indexmz=(r-radius(t)>=vols(i1:i2,j1:j2,k1:k2)); % all points farther
    % distance from closest atom
    vols(i1:i2,j1:j2,k1:k2)...
        =min(vols(i1:i2,j1:j2,k1:k2),r-radius(t));
    % index of closest atom
    vola(i1:i2,j1:j2,k1:k2)...
        =vola(i1:i2,j1:j2,k1:k2).*indexmz+indexm.*t;
    % radius of closest atom
    volr(i1:i2,j1:j2,k1:k2)...
        =volr(i1:i2,j1:j2,k1:k2).*indexmz+indexm.*(radius(t));
    
    % Calculate the EM map caused by protein
    rv=Radius3(2*nw,cfrac1(:,ia));
    volm1=weight(t)*kn*exp(-ka*rv.^2);
    proteinMap(i-nw+1:i+nw,j-nw+1:j+nw,k-nw+1:k+nw)...
        =proteinMap(i-nw+1:i+nw,j-nw+1:j+nw,k-nw+1:k+nw)+volm1;
end;

%% set the mask 0 around and 1 for 1.7 angstrom far away from the protein
% and then shrink the solvent grid mask
disp('shrink');
vols=GridShrink(vols,n,volr,gridfactor);

%% Calculate the solvent density around protein and show the projections and sections
% SolDens=showsection_Func(vola,vols,WaterDensity,n,volm);

% Solvent function model for hydrophilic and hydrophobic surface;
disp('density functions');
if useStepFunction
    dr=.07;
    SolDens=WaterDensity*(1+erf(vols-dr))/2;
else
% Apply different function model to hydrophilic and hydrophobic parts.
    volsC=vols;
    volsC(vola~='C')=-2.65;  % nonpolar
    vols(vola=='C')=-2.65;   % polar
    SolDensE=WaterDensity*((1+erf(vols-0.5))./2+0.2.*exp(-((vols-1.7)./2.5).^2)-0.15.*exp(-((vols-3)./1.5).^2));
    SolDensC=WaterDensity*((1+erf(volsC-1))./2+0.15.*exp(-((volsC-2.2)./2.5).^2)-0.12.*exp(-((volsC-3.6)./1.2).^2));
    SolDens=SolDensE+SolDensC;
end;
compositeMap=proteinMap+SolDens;
disp('done.');
end


%% Generate a contiguous solvent grid mask
function [voln]=GridShrink(vols,n,gridr,gridfactor)

voln=vols;
scaler=1.7;  % Probe radius be used to set the mask 0 around protein within this value.
% indexc=find(vols<=scaler);
voln(vols<=scaler)=-2.65;   %-2.65 is the value which makes EM value caused by solvent smaller than 0.0001.

[X,Y,Z]=ndgrid(1:n,1:n,1:n);

indexgs=find(vols>scaler&vols<=scaler+1.8/gridfactor); % Find a contiguous grid shell around protein;

X1=X(indexgs);Y1=Y(indexgs);Z1=Z(indexgs); % Coordinates of grid shell

cube=Radius3(13,[7 7 7]);  % Core for shrinking
nw2=6;   % Core size around grid shell point, which is used for shrinking

% Shrink from the solvent shell toward protein
for ia=1:length(X1);
    i=X1(ia);
    j=Y1(ia);
    k=Z1(ia);
    
    indexm=cube>gridr(indexgs(ia))+vols(indexgs(ia));
    indexmz=cube<=gridr(indexgs(ia))+vols(indexgs(ia));
    
    voln(i-nw2:i+nw2,j-nw2:j+nw2,k-nw2:k+nw2)...
        =voln(i-nw2:i+nw2,j-nw2:j+nw2,k-nw2:k+nw2).*indexm+indexmz.*vols(i-nw2:i+nw2,j-nw2:j+nw2,k-nw2:k+nw2);
end
end


function r=Radius3(n,org)
% function r=Radius3(n,org)
% Create an n x n x n array where the value at (x,y,z) is the distance from the
% origin.
% The returned values are single precision.
x=org(1);
y=org(2);
z=org(3);
rx=single(1-x:n-x);
ry=single(1-y:n-y);
rz=single(1-z:n-z);
% [X,Y,Z]=ndgrid(1-x:n-x,1-y:n-y,1-z:n-z);  % Make zero at x,y,z
[X,Y,Z]=ndgrid(rx,ry,rz);  % Make zero at x,y,z
r=sqrt(X.^2+Y.^2+Z.^2);
end