% krLookAtMovies

outputDir='jpeg/';
CheckAndMakeDir(outputDir,1);
outN=1920;
outFilt=.1;
clipFraction=1e-4;
figure(1);

d=dir;
nd=numel(d);
for i=3:nd
    name=d(i).name;
    [pa,nm,ex]=fileparts(name);
    if strcmp(ex,'.tif')
        
        disp(['Reading ' name '...']);
        [m0,s]=ReadMovie(name);
        disp('done.');
        ms=RemoveOutliers(sum(single(m0),3));
        %%
        n=NextNiceNumber(size(ms,2));
        ms=Crop(ms,n,0,mean(ms(:)));
        mf=GaussFilt(Downsample(ms,outN),outFilt);
        
        outName=[outputDir nm '.jpg'];
        disp(outName);
        WriteJpeg(mf,outName,clipFraction);
        
        imags(mf);
        title([num2str(i) ':  ' name],'interpreter','none');
        drawnow;
    else
        disp(['Skipping ' name]);
    end;
end;

