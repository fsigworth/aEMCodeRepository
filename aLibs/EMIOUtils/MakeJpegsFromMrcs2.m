% MakeJpegsFromMrcs
% For each *.mrc file in the current directory, make a jpeg and place it
% into the desired directory.
figure(1);
clf;
disp(['Converting mrc images to jpgs in directory ' pwd]);
ds=1;
d=dir;
jpegPath='jpeg_s/';
CheckAndMakeDir(jpegPath,1);

for i=3:numel(d)
    [path, basename, ext]=fileparts(d(i).name);
    if ~d(i).isdir && strcmpi(ext,'.mrc')
        disp(d(i).name);
        [m, s]=ReadMRC(d(i).name);
        pixelsize=s.rez/s.nx;
        imSize=size(m);
            padSize=NextNiceNumber(imSize);
            m=Crop(m-median(m(:)),padSize);
%         m=RemoveOutliers(single(m),4);
        n=size(m,1);
        if ds>1
%             m=BinImage(m,decimation);
            m=Downsample(m,n/ds);
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
