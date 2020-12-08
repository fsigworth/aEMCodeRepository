% MakeJpegsFromMovies

movieDir='/gpfs/gibbs/pi/cryoem/GlaciosYSB/20201203'
outDir='~/scratchs60/20201203/Jpegs/'
fc=.2;
ds=8;
ext='.tif';

CheckAndMakeDir(outDir,1);

d=dir(movieDir);
for i=1:numel(d)
    [pa,nm,ex]=fileparts(d(i).name);
    if strcmp(ex,ext)
        mv=ReadMovie([AddSlash(d(i).folder) d(i).name]);
        disp([d(i).name '  ' num2str(size(mv,3)) ' frames']);
        ms=sum(single(mv),3);
        md=RemoveOutliers(ms);
        % pad and downsample
        mnv=mean(md(:));
        szx=NextNiceNumber(size(md),5,8);
        mdx=Crop(md,szx,0,mnv);
        mds=Downsample(mdx,szx/ds);
        outName=[outDir nm '.jpg'];
        mdsc=imscale(GaussFilt(mds,fc),1,.0001);
        imags(mdsc);
        title(nm,'interpreter','none');
        drawnow;
        WriteJpeg(mdsc,outName);
    end;
end;