%miTrimLog.m
% trim log entries up to the desired point
keepIndex=2; % skipe log entries starting with this
          % e.g. 3: MergeImages
load allNames
mi0=meCreateMicrographInfoStruct14;
nm=numel(allNames);
for i=1:nm
    mi=ReadMiFile(allNames{i});
    seq=miDecodeLog(mi);
    ind=seq(keepIndex);
    if ind>0
        mi.log(seq(keepIndex)+1:end)=[];
    end;
    % Also delete vesicle and partidle information
    mi.vesicle=mi0.vesicle;
    mi.particle=mi0.particle;

    WriteMiFile(mi,allNames{i});
    disp([num2str(seq(cutIndex)) '  ' allNames{i}]);
end;
