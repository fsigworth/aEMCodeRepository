% HRPickingProjs.m
% Make the big library of 2d projections
% Stored directly and as an expansion of eigenimages, with coeffs and
% squared norms of the orginal projections.
%

% cd /Volumes/EMWork/Yangyu/20210224_YW_sel/
% cd ~/EMWork/20210224_YW_sel/
cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Yangyu/20210224_YW_sel') % mini2
% starDir='Refine3D/job110/';
%
refDir='HRPicking/';
ds=3;
zShift=0;

% To use our alpha-subunit composite map
refName='compMap.mrc'
% - Normal:
outTypeString='Comp_'

% - TM centered:
outTypeString='Comp_TMCtr_'
zShift=21;
n=56;

% % To use the TM only reference
% refName='tmMap.mrc';
% n=48;
% zsh=[0 0 21]; % shift TM region to center. This brings it to COM

outDir=[refDir 'Eigs/'];
CheckAndMakeDir(outDir);
eigsBasename=['Eigs' outTypeString];
projsBasename=['projs' outTypeString];

writeProjections=1;
writeEigs=1;

angleSteps=4; % Figure 3 A steps at a radius of 40 A
[map3,s]=ReadMRC([refDir refName]);
pixA=s.pixA*ds;
if zShift ~=0
    sh=[0 0 zShift];
else
    sh=-CenterOfMass(map3)
    end;
    
    map=DownsampleGeneral(circshift(map3,round(sh)),n,1/ds);
ShowSections(map);
%%
angs=hrMakeProjectionAngles(angleSteps,4);
np=size(angs,1);
dotCount=1000;
disp(['Making ' num2str(np) ' projections, ' num2str(dotCount) ' per dot.']);
tic;
projs=rlMakeTemplates(angs,map,1000);
toc;
if writeProjections
    projsFilename=[outDir projsBasename num2str(n) '.mat'];
    disp(['Writing ' projsFilename])
    save(projsFilename,'pixA','map','angs','projs');
end;


%%
disp('Making SVD')
tic
[U,S,V]=svd(reshape(projs,[n*n,np])','econ');
toc;

%
minFracVar=.95;
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

if writeEigs
    CheckAndMakeDir(outDir,1);
    outName=[outDir eigsBasename num2str(n)];
    disp(['Writing ' outName]);
    save(outName,'eigImgs','normCoeffs','norms2','pixA','ds','angs','minFracVar')
    disp('done.');
else
    disp('Eigs not written.');
end;
