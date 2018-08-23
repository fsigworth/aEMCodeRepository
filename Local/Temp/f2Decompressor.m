% f2Decompressor.m
[names, pathName]=uigetfile('*.z.tif','Select compressed files','multiselect','on');
if isnumeric(pathName) % File selection was cancelled
    return
end;
if ~iscell(names)
    names={names};
end;
cd(pathName);
nFiles=numel(names);
for fIndex=1:nFiles
    inName=names{fIndex};
    outName=[inName(1:end-6) '.mrc'];
    disp(['Reading ' inName]);
    disp(['Writing ' outName]);
    
    t=Tiff(inName,'r');
    crf=t.read;
    try
        idString=t.getTag('ImageDescription');
    catch
        error('Not a compressed file');
    end;
    eval(idString);
    s
    close(t);
    toc;
    disp('Decompress');
    tic
    rf=single(crf-128);
    [nx,ny]=size(rf);
    fsModel=ccdAliasedGaussians(nx,s.fs,s.as);
    rm=real(ifftn(fftn(rf).*ifftshift(sqrt(fsModel))));
    ok=1;
    if exist(outName,'file')
        ok=strcmpi(input('File exists.  Overwrite? [n] ','s'),'y');
    end;
    if ok
        WriteMRC(rm,s.pixA,outName);
        disp(' written.');
    else
        disp(' not written.');
    end;
        figure(2);
    subplot(222); imags(GaussFilt(rm,.1));
    subplot(224);
    bf=8;  % bin factor
    semilogy(RadialPowerSpectrum(rm,1,bf));
    drawnow;
end;
% %%
% %
%     figure(2);
%     subplot(221); imags(GaussFilt(md1,.1));
%     subplot(222); imags(GaussFilt(rmr,.1));
%     subplot(224);
%     bf=8;  % bin factor
%     semilogy([RadialPowerSpectrum(md1,1,bf) RadialPowerSpectrum(md1-rm,1,bf) RadialPowerSpectrum(md1-rmr,1,bf)]);
%     %%
%     % subplot(223);
%     % plot([sect(mf) sect(round(mf)-mf)]);
%     subplot(223);
%     plot([sect(rm) sect(md1-rm)]);
% end;


