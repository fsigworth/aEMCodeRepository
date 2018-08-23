% meWriteMergedJpegs
% From a set of mi files, load the corresponding merged images and create
% jpeg copies of them.  Store them in a Merged/Jpeg folder.

doInverseCTF=0;
lfAmp=1;

subtractVesicles=1;
goodVesiclesOnly=0;

ds=1;  % downsampling from existing merged image

[fname, pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
if isnumeric(rootPath)  % user clicked Cancel
    return
end;
cd(rootPath);
%%
figure(1);
SetGrayscale;

nfi=numel(fname);
for fileIndex=1:nfi
    name=fname{fileIndex};
    load([infoPath name]);  % get the mi file
    
    m=meReadMergedImage(mi);
    n=size(m)/ds;
    if ds>1
        m=Downsample(m,n);
    end;
    if subtractVesicles
        if goodVesiclesOnly
            ves=all(mi.vesicle.ok(:,2:3),2);  % vesicle in range and refined
        else
            ves=mi.vesicle.ok(:,1);
        end;
        m=m-meMakeModelVesicles(mi,n,find(ves));
        suffix='mv.jpg';
    else
        suffix='m.jpg';
    end;
    
    if doInverseCTF
        m=meCTFInverseFilter(m,mi,lfAmp);
        dirString='JpegIF/';
        imacs(m);
    else
        dirString='Jpeg/';
        imacs(m);
    end;
    title(mi.baseFilename,'interpreter','none');
    drawnow;
    
    jpegPath=[mi.procPath dirString];
    if ~exist(jpegPath,'dir')
        mkdir(jpegPath);
    end;
    jpegName=[jpegPath mi.baseFilename suffix];
    
    disp(jpegName);
    WriteJpeg(m,jpegName);
end;
