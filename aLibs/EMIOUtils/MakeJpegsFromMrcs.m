% MakeJpegsFromMrcs
figure(1);
clf; SetGrayscale;
disp(['Converting mrc images to jpgs in directory ' pwd]);
ds=1;
d=dir;
jpegPath='';
CheckAndMakeDir(jpegPath,1);

for i=3:numel(d)
    [path, basename, ext]=fileparts(d(i).name);
    if strcmpi(ext,'.mrc')
        disp(d(i).name);
        [m s]=ReadMRC(d(i).name);
        pixelsize=s.rez/s.nx
        imSize=size(m)
        m=RemoveOutliers(single(m),4);
        n=size(m,1);
        if ds>1
%             m=BinImage(m,decimation);
            m=Downsample(m,size(m,1)/ds);
        end;
        WriteJpeg(m,[jpegPath basename '.jpg']);
%         ms=uint8(imscale(m));
%          image(ms);
    imags(m);
        title(d(i).name,'interpreter','none');
        drawnow;
%         imwrite(ms,[jpegPath basename '.jpg'],'jpg');
    end;
end;
disp('Done.');
