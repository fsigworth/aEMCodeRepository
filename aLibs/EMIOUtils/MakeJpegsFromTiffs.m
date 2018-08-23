function MakeJpegsFromTiffs(binning)
% function MakeJpegsFromTiffs(binning)
% For each .tif file in the current directory, make a .jpg file.  The
% optional argument binning controls the degree of downsampling.  Its
% default value is 2.
% fs July 2009

disp(['Converting tiffs to jpgs in directory ' pwd]);
if nargin<1
    binning=2
end;
d=dir;
for i=3:numel(d)
    p=strfind(d(i).name,'.tif');
    if numel(p)>0
        basename=d(i).name(1:p-1);
        disp(d(i).name);
        m=single(ReadEMFile(d(i).name));
        m=RemoveOutliers(m,4);
        n=size(m,1);
        if binning>1
            m=Downsample(m,n/binning);
        end;
%         WriteMRC(single(m), 5/6.7 ,['mrc/' basename '.mrc']);
        ms=uint8(imscale(m));
        imwrite(ms,[basename '.jpg'],'jpg');
    end;
end;
disp('Done.');
