function gFvs=reGroupInsert(ri,iter,iTwin,iGroup,roi,logs)
% function gFvs=reGroupInsert(ri,iter,iTwin,iGroup,roi,logs)
% Returns the Fourier volumes from a group.
useParFor=(ri.nGroups==1);

if iGroup>1
    roiName=[reGetNameCode(ri,iter,iTwin) 'roi.mat'];
    mdisp(logs,[datestr(now) ' Loading the roi file ' roiName]);
    [roi,ok]=reCheckAndLoadMat(ri.outPath, roiName,ri.timeout(2));
    if ~ok
        mdisp(logs,[datestr(now) ' Couldn''nt find the file ' ri.outPath roiName]);
        error(['Couldn''nt find the file ' ri.outPath roiName]);
    end;
%     Quick check whether we've really loaded the whole thing
    if ~isfield(roi,'classNorms') || ~isfield(roi,'classMeans') % didn't get the whole struct
        mdisp(logs,[datestr(now) ' Trying again to load the roi struct ']);
        pause(5);
        [roi,ok]=reCheckAndLoadMat(ri.outPath, roiName,10);
        if ~ok
            mdisp(logs,[datestr(now) ' Couldn''nt load the entire roi struct ' ri.tempPath roiName]);
            error(['Couldn''nt load the entire roi struct ' ri.tempPath roiName]);
        end;
    end;
end;

nRefs=size(ri.angles,1);
rk=struct;  % make a fake ri structure for reGetGroupFlags
rk.nTwin=nRefs;
rk.nGroups=ri.nGroups;

% Select the class means and norms for this group
gFlags=reGetGroupFlags(rk,1,iGroup);  % get the flags for splitting up the classes
ngRefs=sum(gFlags);
classNorms=roi.classNorms(:,:,gFlags,:);
classMeans=roi.classMeans(:,:,gFlags,:);
angles=ri.angles(gFlags,:);

% % %%%%%%
% % classNorms=Downsample(classNorms,48,1);
% % classMeans=Downsample(classNorms,48,1);

mdisp(logs,[datestr(now) ' Reconstruction start']);
n=size(classNorms,1);
gProjs=zeros([n n ngRefs 2 ri.nVols],'single');  % The big array to give to Fourier Insertion

for iVol=1:ri.nVols
    realNorms=zeros(n,n,ngRefs,'single');
    for i=1:ngRefs
        realNorms(:,:,i)=fftshift(real(ifftn(ifftshift(classNorms(:,:,i,iVol)))));
    end;
    gProjs(:,:,:,1,iVol)=classMeans(:,:,:,iVol);
    gProjs(:,:,:,2,iVol)=realNorms;
end;
gFvs=reFourierInsertion2(gProjs,angles,ri.symmetry,useParFor);  % no parfor

% Save the group Fourier volumes
fvName=[reGetNameCode(ri,iter,iTwin,iGroup) 'fvs.mat'];
mdisp(logs,[datestr(now) ' saving Fourier vols: ' fvName]);
% save([ri.tempPath fvName],'gFvs');
save([ri.tempPath fvName '_tmp'],'gFvs');
eval(['!mv ' ri.tempPath fvName '_tmp ' ri.tempPath fvName]);

