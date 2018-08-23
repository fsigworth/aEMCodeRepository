function [gRoi,refs]=reNodeEMStep(ri,gSi,gMoi,gImgs,gAltImgs,logs)
% Given the group si, moi and images, perform an EM step and return the
% run-output structure roi and the class means cls.  Internally, we split
% the job into ri.nSlices slices and run them in parallel, collecting the
% results before returning.

% Make the flags for the slices
nGroup=size(gImgs,3);

sliceFlags=reGetSliceFlags(ri,nGroup);
ns=ri.nSlices;

ssi=cell(ns,1);
simgs=cell(ns,1);
saltImgs=cell(ns,1);
sRoi=cell(ns,1);

for i=1:ns
    [ssi{i},simgs{i},saltImgs{i}]=rsStackSplit(sliceFlags(:,i),gSi,gImgs,gAltImgs);
end;

% disp('Making refs');
if ri.flags.mode2D % 2D classification
    refs=gMoi.refs;
else
    refAngles=reGetAngleList(ri,false);  % pick beta, gamma angles
    refs=reMakeTemplates(gMoi.refVols,refAngles);  % refs(x,y,iRef,iVol)
end;
mdisp(logs,['  ' num2str(size(refs,3)) ' references']);

if ~isfield(gMoi,'activeAlphas')
    gMoi.activeAlphas=true([size(ri.alphas) size(gImgs,3)]);
end;
% ------------------------ Do the EM step ---------------------------
% ------------------------------------------------------------------
    for is=1:ns  % Loop over slices
        if ns>1
            disp(is);
        end;
        mo1=reSplitStructureFields(gMoi,sliceFlags(:,is));
        sRoi{is}=reEMStep26(simgs{is},refs,ssi{is},ri,mo1,saltImgs{is},logs);
        if ri.flags.removeRings
            mRefs=reMakeMatchedRefs(refs,ssi{is},ri,sRoi{is});
            sRoi{is}.matchedRefs=mRefs;
            sRoi{is}.ringFits=reFitVesicleRings(ssi{is},simgs{is}-mRefs,...
                              ri.ringFitMask,ri.ringRadii,ri.ringWidths);
        else
            sRoi{is}.ringFits=0;
        end;
    end;
% ------------------------------------------------------------------
% ------------------------------------------------------------------

%% -----Accumulate slices----------
gRoi=reCollectStructureFields(sRoi,sliceFlags);

