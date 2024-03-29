% HRPickingProjs.m
% Make the big library of 2d projections
% Stored directly and as an expansion of eigenimages, with coeffs and
% squared norms of the orginal projections.
%

% HR picking tests
% cd /Volumes/EMWork/Yangyu/20210224_YW_sel/
% cd ~/EMWork/20210224_YW_sel/
cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Yangyu/20210224_YW_sel') % mini2
refDir='HRPicking/';

% Kv dataset
cd('/Volumes/D255/20181216/No5Graphene/sq05_1/');
inDir='HRRefs/';
outDir='HRRefs/';
%
ds=3;
zShift=0;
centroidShift=1; % if zShift is zero, move to centroid.

% % To use our alpha-subunit composite map
% refName='compMap.mrc'
% outTypeString='Comp_';
% n=56;
% zShift=-10;

% Weak TM domain
refName='compWeakMap.mrc'
outTypeString='CompWeakUnSh_';
n=56;
centroidShift=0;

% % - TM centered:
% outTypeString='Comp_TMCtr_'
% zShift=21;
% n=56;

% % To use the TM only reference
% refName='tmMap.mrc';
% n=48;
% outTypeString='TM_';
% zsh=[0 0 21]; % shift TM region to center. This brings it to COM

angleSteps=[4 4 4];
% angleSteps=[2.5 2.5 4];
stepString=num2str(angleSteps(1));
for i=2:3
    stepString=[stepString '_' num2str(angleSteps(i))];
end;
disp(stepString);


CheckAndMakeDir(outDir);
eigsBasename=['Eigs' outTypeString];
projsBasename=['projs' outTypeString];

writeProjections=1;
makeEigs=0;
writeEigs=0;

% angleSteps=4; % --produces 3 A steps at a radius of 40 A
[map3,s]=ReadMRC([inDir refName]);
pixA=s.pixA*ds;
if zShift ~=0
    sh=[0 0 zShift];
elseif centroidShift
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
    projsFilename=[outDir projsBasename '_' stepString '_' num2str(n) '.mat'];
    disp(['Writing ' projsFilename])
    save(projsFilename,'pixA','map','angs','projs');
end;


%%
if makeEigs
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
end;