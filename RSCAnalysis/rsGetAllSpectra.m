% rsGetAllSpectra.m
movieMode=1;

% Have the user select some mi files
[fname, pa]=uigetfile({'*mi.txt' '*mi.mat'},'Select mi files to examine','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;

%%
ds=8;
% pick up all the mi structures and compute spectra
[basePath,infoPath]=ParsePath(pa);
cd(basePath);
nim=numel(fname);


%%
for i=1:nim
    
    disp(fname{i});
    mi=ReadMiFile([infoPath fname{i}]);
    mis{i,1}=mi;
    if i==1
        n=mi.imageSize;
        ssp=zeros(n);
        spMask=1-fuzzymask(n/ds,2,n(1)/(8*ds),0);
    end;
    
    if movieMode
        frame1=mi.frameSets(1,1);
        frameN=mi.frameSets(1,2)-frame1+1;
        movieName=[mi.moviePath mi.movieFilename];
        if exist(movieName,'file')
            disp(['  ' movieName]);
            mvd=double(ReadMovie(movieName,frame1,frameN));
            mvd=min(mvd,10);  % restrict to small counts
            mvMean=mean(mvd(:));
            mvd=Crop(mvd,mi.imageSize,1,mvMean);
            
            for j=1:frameN
                sp=fftshift(abs(fftn(mvd(:,:,j))).^2);
                ssp=ssp+sp;
            end;
            bsp=BinImage(sqrt(ssp),ds);
            msp=spMask(:)'*bsp(:)/(spMask(:)'*spMask(:));
            bsp=spMask.*bsp+(1-spMask)*msp;
            imags(bsp);
            title(i);
            drawnow;
            
        else
            disp(' --not found');
        end;
    else
        iName0=[mi.imagePath mi.imageFilenames{1}];
        [imageName,ok]=CheckForImageOrZTiff(iName0);
        if ok
            disp(['  ' imageName]);
            [m,pixA,ok]=ReadEMFile(imageName);
            if ok
                sp=fftshift(abs(fftn(double(m))).^2);
                ssp=ssp+sp;
                bsp=BinImage(sqrt(ssp),ds);
                msp=spMask(:)'*bsp(:)/(spMask(:)'*spMask(:));
                bsp=spMask.*bsp+(1-spMask)*msp;
                imags(bsp);
                title(i);
                drawnow;
            end;
        else
            disp('  skipped');
        end;
    end;
end;
%%
if movieMode
    save('MovieSpectrum.mat','ssp','bsp','pixA');
else
    save([mi.imagePath 'SpectraSum.mat'],'ssp','bsp','pixA');
end;