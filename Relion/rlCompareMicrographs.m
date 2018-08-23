% rlCompareMicrographs.m
% Check to see that two sets of micrographs agree.

ncc=64; % cropped cc size

[miNames, pa]=uigetfile({'*mi.txt' '*mi.mat'},'Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    cd(rootPath);
if ~iscell(miNames)
    miNames={miNames};
end;
%%
disp('Getting micrograph path 1');
miPath1=AddSlash(uigetdir);

disp('Getting micrograph path 2');
miPath2=AddSlash(uigetdir);

figure(1);

nim=numel(miNames);
step=100;
%%
for i=1:step:nim
    mi=ReadMiFile([infoPath miNames{i}]);
    for k=1:numel(mi.imageFilenames)
    [nm1,ok1]=CheckForImageOrZTiff([miPath1 mi.imageFilenames{k}]);
    [nm2,ok2]=CheckForImageOrZTiff([miPath2 mi.imageFilenames{k}]);
    if ok1 && ok2
        m1=ReadEMFile(nm1);
        m2=ReadEMFile(nm2);
        cc=Crop(fftshift(real(fftn(m1).*conj(fftn(m2)))),ncc);
        rms1=std(m1(:));
        rmsd=std(m1(:)-m2(:));
        mb1=BinImage(m1,4);
        mb2=BinImage(m2,4);
        [ms1,mulr,addr]=imscale(mb1,256,.001);
        subplot(221);
        imaga(ms1);
                title(nm1,'interpreter','none');
        subplot(222);
        imaga(mb2*mulr+addr);
                title(nm2,'interpreter','none');
        subplot(223);
        plot([sect(mb1) sect(mb2) sect(mb1-mb2)]);
        subplot(224);
        imags(cc);
        title(num2str([rmsd rms1]));
        drawnow;
    else
        disp([num2str([ok1 ok2]) ' missing: ' mi.imageFilenames{k}]);
    end;
    end;
end;
