function roi=reGatherRois(ri,iter,iTwin,logs,writeRoi)
% function roi=reGatherRois(ri,iter,iTwin,logs,writeRoi)
% Read the group roi files and combine them, returning the composite roi.
% function  ok=reGatherRois(ri,iter,iTwin)
% Just check for the existence of the files

if nargin<4 % we just check for existence of the files
    roi=true;  % output varible will be boolean
    for ig=1:ri.nGroups
        roiName=[reGetNameCode(ri,iter,iTwin,ig) 'roi.mat'];
        roi=roi && exist([ri.tempPath roiName],'file');
    end;
    return;
    
else
    if nargin<6
        writeRoi=true;
    end;
        groupFlags=reGetGroupFlags(ri,iTwin);
        rois=cell(ri.nGroups,1);
%         rois{1}=gRoi;
        mdisp(logs,[datestr(now) ' Gathering roi files']);
        for ig=1:ri.nGroups
            roiName=[reGetNameCode(ri,iter,iTwin,ig) 'roi.mat'];
            mdisp(logs,[datestr(now) ' reading group: ' roiName]);
            [gRoi,ok]=reCheckAndLoadMat(ri.tempPath, roiName,ri.timeout(2),1,1);
            if ~ok
                mdisp(logs,['Couldn''t find the file ' ri.tempPath roiName datestr(now)]);
                error(['Couldn''nt find the file ' ri.tempPath roiName]);
            end;
            rois{ig}=gRoi;
        end;
        mdisp(logs,[datestr(now) '  all found.']);
        roi=reCollectStructureFields(rois,groupFlags);
    if writeRoi
        roiName=[ri.outPath reGetNameCode(ri,iter,iTwin) 'roi.mat'];
        mdisp(logs,[datestr(now) ' Saving the roi file ' roiName]);
        save([roiName '_tmp'],'roi');
        eval(['!mv ' roiName '_tmp ' roiName]);
    end;
    
end;
