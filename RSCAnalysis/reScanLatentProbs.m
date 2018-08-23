% reScanLatentProbs.m
% Read the *moi.mat files from RSC-ML reconstruction and display the latent
% probabilities.

d=dir;
ir=0;
im=0;
as=zeros(0,3);
for ind=3:numel(d)
        a=sscanf(d(ind).name,'%c%d%s'); % e.g. i01a_moi.mat -> i 1 a_moi.mat
        if numel(a)==11 && a(5)=='m' % got all fields
            im=im+1;
            ams(i,:)=[a(2) a(3)-96 ind];
        elseif numel(a)==1 && a(5)=='r' % roi
            ir=ir+1;
            ars(i,:)=[a(2) a(3)-96 ind];
        else
            a2=sscanf(d(ind).name,'%c%c
            return
        end;
end;
firstIter=min(as(:,1))
lastIter=max(as(:,1))
load(d(as(1,3)).name);
pVols=moi.pVols';
nVols=numel(pVols);
iters=as(1,1);
if iters==1
    pVols(1,:)=NaN;   % don't plot the first point.
end;
for i=2:size(as,1)
    name=d(as(i,3)).name;
    disp(name);
    load(name);
    pVols(i,:)=moi.pVols';
    iters(i,1)=as(i,1);
end;
%%



if iter>ri.startIter && iGroup>1  % read moi from a file
        mdisp(logs,[datestr(now) ' Waiting for the moi file ']);
        [moi,ok]=reCheckAndLoadMat(ri.outPath, [reGetNameCode(ri,iter,iTwin) 'moi.mat'],...
            ri.timeout(1),1,numel(fieldnames(moi)));


        if iGroup>1 || ri.flags.writeAllGroupFiles  % Encode iter, iGroup and iTwin in the output file name
            roiName=[ri.tempPath reGetNameCode(ri,iter,iTwin,iGroup) 'roi.mat'];
            save([roiName '_tmp'],'gRoi');
            eval(['!mv ' roiName '_tmp ' roiName]);
        end;

        %%
        if showGraphics
            %     Display the group's latent probabilities
            figure(1);
            set(gcf,'position',[100,100,400,360]);
            figure(2);
            set(gcf,'position',[500,100,800,720]);
            reShowLatentVars(imgs,refs,ri,gMoi,gRoi,iter,[1 2]);
            drawnow;
            if ri.flags.saveFigures
                set(gcf,'paperpositionmode','auto');
                figName=[reGetNameCode(ri,iter,iTwin,iGroup) '.jpg'];
                print([ri.outPath 'jpeg/' figName],'-djpeg');
            end;
        end;
end;