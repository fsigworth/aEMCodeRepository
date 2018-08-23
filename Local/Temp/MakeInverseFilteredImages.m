% MakeInverseFilteredImages


    [names, pathName]=uigetfile('*mi.*','Select mi files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(names)
        names={names};
    end;
    [rootPath, infoPath]=ParsePath(pathName);
    cd(rootPath);
    miNames=cell(0);
    for i=1:numel(names)
        miNames{i}=[infoPath names{i}];
    end;
%%
    
    lfAmp=.5 ;
    fDeTrend=.00002;
    fHP=.00002 ;
    fLP=.2; 
    
    figure(1); SetGrayscale;
    nim=numel(miNames);
    for i=1:nim
        mi=ReadMiFile(miNames{i});
        m=ReadMRC([mi.procPath mi.baseFilename 'm.mrc']);
        m1=Downsample(m,size(m)/2);
            mfilt=meCTFInverseFilter(m1,mi,lfAmp,fDeTrend,fHP);
        imacs(GaussFilt(mfilt,fLP));
        title(miNames{i},'interpreter','none');
         pause;
    end;
    disp('Done.');
    