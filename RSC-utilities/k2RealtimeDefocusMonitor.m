% k2RealtimeDefocusMonitor

nmiMin=1;
nMiHistory=1000;
maxJobs=100;
doRsync=0;

rawDir='/ysm-gpfs/pi/cryoem/krios/20170808/20170811/';
expDir='~/scratch60/170808/KvLipo/';
movieDir='movie_frames/sq03/';
ok=true;
while ok
if doRsync
    rsyncStr=['rsync -tv ' rawDir '* ' expDir movieDir];
    disp(rsyncStr);
    s=system(rsyncStr)
    
    cd(expDir);
    k2CreateInfoFiles;
    nmi=numel(miNames);
    disp([num2str(nmi) ' new files.']);
%%
    if nmi>= nmiMin
        njobs=min(nmi,maxJobs);
        allNames=miNames;
        save('allNames.mat','allNames');
        cd ~/
        execStr=['sbatch --array=1-' num2str(njobs) ' k2pRun'];
        disp(execStr);
        system(execStr);
    end;
end;
%%    
    cd(expDir);
    d=dir('Info/');
    nMis=0;
    miNames=cell(1,0);
    for i=1:numel(d)
        nm=d(i).name;
        q=strfind(nm,'mi.txt');
        if numel(q)==1 % valid mi file
            nMis=nMis+1;
            miNames{nMis,1}=['Info/' nm];
        end;
    end;
%%
if nMis>0
        miStart=max(nMis-nMiHistory+1,1);
        miNames(1:miStart-1)=[];
    figure(1);
    
    nMis=numel(miNames);
    mis=cell(nMis,1);
disp(['Reading ' num2str(nMis) ' mi files']);
    for i=1:nMis
%        disp(miNames{i});
        mis{i}=ReadMiFile(miNames{i});
    end;
    disp('done');
    miDefocusTrend;
    grid on;
    end;
    return
       pause(30);     
end;

