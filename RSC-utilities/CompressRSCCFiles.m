% CompressRSCCFiles
clear;
qN=500;

[fname, pa]=uigetfile('*mi.*','Select mi files','multiselect','on');
if isnumeric(pa) % File selection cancelled
    return
end;
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(rootPath);

for i=1:numel(fname)
    mi=ReadMiFile([infoPath fname{i}]);
    rsName=[mi.procPath mi.baseFilename 'rscc.mat'];
    if exist(rsName,'file')
        load(rsName);
        mVesGood=round(qN*mVesGood)/qN;
        mVesBad=round(qN*mVesBad)/qN;
        mxCC=round(qN*mxCC)/qN;
        mxVars=round(mxVars);
        mxDist=round(mxDist);
        vList=single(vList);
        angleList=single(angleList);
        clear gMaxVals;
        if ~exist('log','var')
            log='';
        end;
        save(rsName,'mxCC','mxVars','mxVesInds',...
            'mxDist','mxTemplInds','mxRsos','partRadius', 'membraneOffsetA','ds',...
            'badVesMask','eigenImgs','vList','angleList','ppVals','pwfRef',...
            'mVesGood','mVesBad','log');
        disp(['written:   ' mi.procPath mi.baseFilename 'rscc.mat']);
    else
        disp(['not found: ' mi.procPath mi.baseFilename 'rscc.mat']);
    end;
end;