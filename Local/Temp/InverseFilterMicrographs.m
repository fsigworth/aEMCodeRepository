% InverseFilterMicrographs

[names, pathName]=uigetfile('*mi.*','Select mi files','multiselect','on');
if isnumeric(pathName) % File selection was cancelled
    return
end;
if ~iscell(names)
    names={names};
end;
[rootPath, infoPath]=ParsePath(pathName);
cd(rootPath);
%%
nd=1024;

nFiles=numel(names);
imgsUnfilt=zeros(nd,nd,nFiles,'single');
imgsFilt=zeros(nd,nd,nFiles,'single');

for iFile=1:nFiles
    miName=[infoPath names{iFile}];
    mi=ReadMiFile(miName);
    doses=mi.doses
    %     Get the merged image
    m=ReadMRC([mi.procPath mi.baseFilename 'm.mrc']);
    md=Downsample(m,nd);
    subplot(221);
    imags(GaussFilt(md,.1));
    title(mi.baseFilename,'interpreter','none');
    subplot(223);
    mdo=(md+1)*mi.doses(1);
    plot(mean(mdo,2));
    axis([0 1024 0 45]);
    imgsUnfilt(:,:,iFile)=mdo;
    drawnow;
    %%
    for i=1:numel(mi.ctf)
        mi.ctf(i).alpha=.035;
    end;
    [mf,H]=meCTFInverseFilter(md,mi,1,0,0);
    subplot(222);
    title(iFile);
    imags(mf);
    subplot(224)
    plot(sect(mf*mi.doses(1)*mi.doses(1)));
    drawnow;
    imgsFilt(:,:,iFile)=mf*mi.doses(1)+mi.doses(1);
    %%
end;