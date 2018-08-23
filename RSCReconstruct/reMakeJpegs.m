% reMakeJpegs.m
% Create graphics files (ShowSections of models) in case they weren't
% created by reEMReconstruct3

doAlignVolumes=1;
fscRadius=.22;
fscEdge=.2;
fscZ=0;
startIter=29;
writeJpegs=1;
showMask=0;

% [riName, pa]=uigetfile('*ri.mat','Select the ri.mat file');
riName='ri.mat';
pa=uigetdir('.','Select a Recon directory');

if isnumeric(pa)
    return
end;

cd(pa);
[rootDir,localDir]=ParsePath(pa);
localDir
if ~exist(riName,'file')
    display('No ri.mat file found');
    return
end;
CheckAndMakeDir('jpeg/');
load(riName);
figure(3);
clf;
pos=get(gcf,'position');
set(gcf,'position',[pos(1:2),630,600]);
n=1;
tx=0;
gamma=0;
mirror=0;
disMsk=1;

% Create jpegs from mrc volumes
if writeJpegs
    disp('Writing jpeg files...');
end;
%%
for iter=startIter:ri.nIters
    for iVol=1:ri.nVols
        twinVols=[];
        msk=1;
        twinPresent=zeros(ri.nTwins,1);
        for iTwin=1:ri.nTwins
            volName=['mrc/' reGetNameCode(ri,iter,iTwin,-iVol) '.mrc'];
            if exist(volName,'file')
                [vol,s]=ReadMRC(volName);
                twinVols(:,:,:,iTwin)=vol;
                twinPresent(iTwin)=1;
            end;
        end;
        if sum(twinPresent)>1
            n=size(twinVols,1);
            freqs=(.5:n/2-.5)/(n*s.pixA);
            fscZ=-.13;
            fscRadius=.2;
            msk=fuzzymask(n,3,n*fscRadius,n*fscEdge,ceil((n+1)/2)+[0 0 n*fscZ]);
            fscZ=.13;
            fscRadius=.2;
            msk=min(1,msk+fuzzymask(n,3,n*fscRadius,n*fscEdge,ceil((n+1)/2)+[0 0 n*fscZ]));
            
            if doAlignVolumes
                [twinVols(:,:,:,2),tx,gamma,mirror]=reAlignVolumes(msk.*twinVols(:,:,:,1),twinVols(:,:,:,2));
            end;
            
            fsc=FSCorr2(msk.*twinVols(:,:,:,1),msk.*twinVols(:,:,:,2));            
        else
            fsc=0;
        end;
        if showMask
            disMsk=msk;
        else
            disMsk=fuzzymask(n,3,n*.45,n*.1);
        end;
        if any(twinPresent)
            twinVols(:,:,:,ri.nTwins+1)=sum(twinVols(:,:,:,twinPresent>0),4);
            twinPresent(ri.nTwins+1)=1;
            for iTwin=1:ri.nTwins+1
                if twinPresent(iTwin) || iTwin>ri.nTwins+1
                    ShowSections2(disMsk.*twinVols(:,:,:,iTwin),[],45);
                    subplot(3,3,1);
                    title(localDir);
                    figName=[reGetNameCode(ri,iter,iTwin,-iVol) '.jpg'];
                    subplot(3,3,2);
                    title(figName);
                    if sum(twinPresent)>1
                        subplot(3,3,9);
                        fsc(1)=1;
                        fth=.143;
                        plot(freqs,[fsc 0*fsc+.143]);
                        axis([0 .15 0 1]);
                        p0=find(fsc<fth,1);
                        if numel(p0)<1
                            p0=inf;
                        end;
                        p0=min(p0,numel(fsc)-1);
                        p=p0;
                        if fsc(p0-1)>fsc(p0)
                            p=p0-(fth-fsc(p0))/(fsc(p0-1)-fsc(p0));
                        end;
                        res=(p-.5)/(n*s.pixA);
                        xlabel('Frequency, Å^{-1}');
                        ylabel('FSC');
                        title(num2str([tx gamma mirror res],3));
                    end;
                    drawnow;
                    if writeJpegs
                        set(gcf,'paperpositionmode','auto');
                        print(['jpeg/' figName],'-djpeg');
                    end;
                    %                     disp(figName);
                end;
            end;
        end;
    end;
end;
disp(['last file: ' figName]);
disp('done.');