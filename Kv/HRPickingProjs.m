% HRPickingProjs.m
% Make the big library of 2d projections
% Stored as an expansion of eigenimages, with coeffs and squared norms of
% the orginal projections.
%

% cd ~/EMWork/Yangyu/20210224_YW_sel/
% cd /Volumes/EMWork/Yangyu/20210224_YW_sel/
cd ~/EMWork/20210224_YW_sel/
% starDir='Refine3D/job110/';
%
% To use our alpha-subunit composite map
refDir='HRPicking/';
refName='tmMap.mrc';
zsh=[0 0 21]; % shift TM region to center

outDir=[refDir 'Eigs/'];
outBasename='EigsTM_'
writeProjections=1;
projsFilename=[outDir 'projsTm.mat'];

doWrite=0;
ds=3;
n=48
angleSteps=4; % Figure 3 A steps at a radius of 40 A
[map3,s]=ReadMRC([refDir refName]);
pixA=s.pixA*ds;
map=DownsampleGeneral(circshift(map3,zsh),n,1/ds);
ShowSections(map);
%%
angs=hrMakeProjectionAngles(angleSteps,4);
np=size(angs,1);
disp('Making Projections');
tic;
projs=rlMakeTemplates(angs,map,1000);
toc;
if writeProjections
    save(projsFilename,'pixA','map','angs','projs');
end;


%%
disp('Making SVD')
tic
[U,S,V]=svd(reshape(projs,[n*n,np])','econ');
toc;

%
minFracVar=.98;
nv=n*n;
coeffs=U*S;
allVars=cumsum(coeffs.^2,2);
norms2=allVars(:,end);
vars=allVars./repmat(norms2,1,nv);

minVars=min(vars',[],2);
nc=find(minVars>minFracVar,1); % find the number of terms to bring us within 1%
disp(['nc = ' num2str(nc) ' for ' num2str(100*minFracVar) '% variance']);
% maxVars=max(vars',[],2);
% plot([minVars maxVars]);

coeffs=coeffs(:,1:nc);
normCoeffs=coeffs./repmat(sqrt(norms2),1,nc);
eigVs=V(:,1:nc);
eigImgs=reshape(eigVs,n,n,nc);
recProjs=reshape(eigVs*normCoeffs',[n,n,np]);

figure(3);
for i=1:8,mysubplot(4,4,i); imags(projs(:,:,i*1000)); end;
for i=1:8,mysubplot(4,4,i+8); imags(recProjs(:,:,i*1000)); end;

if doWrite
    CheckAndMakeDir(outDir,1);
    outName=[outDir outBasename num2str(n)];
    disp(['Writing ' outName]);
    save(outName,'eigImgs','normCoeffs','norms2','pixA','ds','angs')
    disp('done.');
else
    disp('Nothing written.');
end;
