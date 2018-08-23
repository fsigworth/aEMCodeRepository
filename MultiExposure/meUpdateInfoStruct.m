function mo=meUpdateInfoStruct(mi,newVersion)
if nargin<2
    newVersion=10;
end;
mo=mi;
mo.version=newVersion;
if mi.version==9 %  We will convert from this
    if ~isfield(mo,'imagePath')
        mo.imagePath=mi.localPath;
    end;
    if ~isfield(mo,'procPath')
        mo.procPath=mi.localPath;
    end;
    if ~isfield(mo,'infoPath')
        mo.infoPath=mi.localPath;
    end;
    if ~isfield(mo,'imageFilenames')
        mo.imageFilenames=mi.rawFile;
    end;
    if ~isfield(mo,'imageA')
        if numel(mi.nPix)<2
            mi.nPix=mi.nPix*[1 1];
        end;
        mo.imageA=single(mi.nPix)*mi.pixA;
    end;
    mo=rmfield(mo,{'pixA','nPix','nExposures','localPath'});
end;
