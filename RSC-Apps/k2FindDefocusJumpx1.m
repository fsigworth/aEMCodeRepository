% k2FindDefocusJump(gInfoPath,gMiNames)

% Attempt to find the frame where the defocus jump occurs in a K2 movie.
% 
% Sets the following fields in the mi file:
% pixA
% frameDose
% frameSets
pars=struct;

pars.pixADefault=1.247;
pars.defoci=[1.5 10];
pars.testSegs=[2 16; 30 inf];
pars.ds=4;
pars.ds1=4;
pars.fo=1/10;
pars.fi=1/200;
pars.doSaveFigures=1;

forceNewMiFrameSets=0;

nFiles=numel(gMiNames);


if nFiles>1
    disp(['Working on ' num2str(nFiles) ' files.']);
end;
nSteps=size(pars.testSegs,1)-1;
gTransFrames=zeros(nFiles,nSteps);

for ind=1:nFiles
% parfor ind=1:nFiles
    infoName=[gInfoPath gMiNames{ind}];
    miStart=load(infoName);
    if forceNewMiFrameSets || numel(miStart.mi.frameSets)<2
        mi=k2FindJumpFcn(miStart.mi,pars);
        fullName=meSaveMiFile(mi);
        %%
        disp(['Updated: ' fullName]);
        gTransFrames(ind,:)=mi.frameSets(1:end-1,2)'+1;
        disp(' ');
    else
        disp(['Skipped:  ' infoName]);
        disp(' ');
    end;
end;
%%
% Show the final results
nr=3;
nc=4;
subplot(nr,nc,nr*nc-1);
hist(gTransFrames);
disp(' ');
disp('Summary:');
for i=1:nFiles
    disp([gMiNames{i} '   ' num2str(gTransFrames(i,:)')]);
end;
