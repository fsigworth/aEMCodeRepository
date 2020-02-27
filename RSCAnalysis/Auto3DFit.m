% Auto3DFit.m

% Get the map to dock
% fc=1;
% inName='ExpMono';
% inMatrix='ExpMonoDocking1_L.mat';
% refScale=4;

fc=.025;
inName='PdbMono';
inMatrix='PdbMonoDocking1_L.mat';
refScale=10;

[rawMap,s]=ReadMRC([inName 'Map.mrc']);
m0=GaussFilt(rawMap,fc*s.pixA);
load(inMatrix); % loads mat
mat0=mat; % manually-determined matrix.

% Get the reference map
[mRef,s]=ReadMRC('DimerMap.mrc');


% load refsmModel.mat % load refs
% me=GaussFilt(refsm,.2);
% load refMatrix.mat


m1=ERotate3(m0,mat0);
mRef=mRef*refScale;
n=size(m1,1);
n2c=ceil((n+1)/2);

% Fit map m1 into m2
r2d=mRef(:);
f2c=conj(fftn(mRef));

P=zeros(1,3);
steps=0.03*ones(1,3);
shifts=zeros(3,1);

p=Simplex('init',P,steps);
figure(1);

% return
% 
nIters=50;
for i=1:nIters
    mat1=[EulerMatrix(P) [0;0;0]; [0 0 0 1]];
    r1=ERotate3(m1,mat1,[],1);
    ccf=real(ifftshift(ifftn(fftn(r1).*f2c)));
    [val,sx,sy,sz]=max3d(ccf);
    r1sh=circshift(r1,n2c-[sx sy sz]);
    if mod(i,20)==0
    disp([NCC P]);
        ShowSections(r1sh+mRef);
        title(NCC);
        drawnow;
    end;
    r1d=r1sh(:);
    NCC=r1d'*r2d/sqrt(r1d'*r1d);
    P=Simplex(-NCC);
end;
%%
% r1sh is now the best transformed map.

ShowSections(r1sh+mRef);
title('Last fit iteration');
drawnow;

%% Final transformation
P=Simplex('centroid');
    ccf=real(ifftshift(ifftn(fftn(r1).*f2c)));
    [val,shifts]=max3di(ccf);
% % % % % % % % %     doesn't work.
    
    mat1Final=[EulerMatrix(P) n2c-shifts'; [0 0 0 1]];
% %     matFit=[EulerMatrix(P) [0;0;0] ; [0 0 0 1]];
%     mapFit=real(ifftn(ERotate3(m1,matFit)).*S);
    compositeMat=mat1Final*mat0;
    mapFit=ERotate3(m0,compositeMat);
    clf;
    ShowSections(mapFit+mRef);
    title('Final matrix');
    drawnow;
    disp(compositeMat);
% %     ShowSections(mapFit+m2);
if mat0(1,4)>0 % shift to right
    side='R';
else
    side='L';
end;
outName=[inName 'Docked_' side '.mrc'];
disp(['Writing ' outName]);
WriteMRC(mapFit,s.pixA,outName);
disp('done.');
