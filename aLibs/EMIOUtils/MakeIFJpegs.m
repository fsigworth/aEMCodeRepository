% MakeIFJpevgs.m
% Read each mi file, and if we find the merged image, do inverse filtering and write
% out a jpeg into the desired directory.
% We assume we're in the experiment directory.

displayOn=1;
targetSize=960;
suffixes={'s' 'vs'};
inverseFraction=.35;
lpFilter=.05;
jpegPath='jpeg/';

if displayOn
    figure(1);
    clf;
end;
disp(['Creating inverse-filtered jpgs in directory ' jpegPath ' in ' pwd]);
CheckAndMakeDir(jpegPath,1);
miNames=f2FindInfoFiles;
% first=1;
for i=1:numel(miNames)
    mi=ReadMiFile(miNames{i});
    for j=1:numel(suffixes)
        [m,name,ok]=meReadMergedImage(mi,0,suffixes{j});
        if ok
            disp(name);
            ds=round(size(m,1)/targetSize);
            approxPixA=mi.pixA*mi.imageSize(1)/outputSize(1);
            outputSize=round(size(m)/ds);
            md=DownsampleGeneral(m,outputSize,1/ds);
            [mf,H]=meCTFInverseFilter(md,mi,inverseFraction);
            mf2=GaussFilt(mf,lpFilter*approxPixA);
            [pa,nm,ex]=fileparts(name);
            suffix=sprintf('.w%02d',100*inverseFraction);
            outName=[jpegPath nm suffix '.jpg'];
            WriteJpeg(mf2,outName);
            if displayOn
                imags(mf2);
                title(outName,'interpreter','none');
                drawnow;
            end;
        end;
    end;
end;
% 
% d=dir;
% jpegPath='../jpeg/';
% CheckAndMakeDir(jpegPath,1);
% extension='.mrc';
% suffix='';
% saveUnfiltered=1;
% saveInverseCTF=.3; % if nonzero, fraction of inverse filtering.
% 
% for i=1:numel(d)
%     [path, basename, ext]=fileparts(d(i).name);
%     if ~d(i).isdir && strcmpi(ext,extension) && strndcmp(basename,suffix,numel(suffix))
%         disp(d(i).name);
%         [m, s]=ReadMRC(d(i).name);
%         pixelsize=s.rez/s.nx;
%         imSize=size(m);
% %         m=RemoveOutliers(single(m),4);
%         n=size(m,1);
%         if ds>1
% %             m=BinImage(m,decimation);
%             m=Downsample(m,size(m,1)/ds);
%         end;
% if saveUnfiltered
%     WriteJpeg(m,[jpegPath basename '.jpg']);
%     imags(m);
%     title(d(i).name,'interpreter','none');
%     drawnow;
% end;
% if saveInverseCTF>0
%     mf=meCTFInverseFilter(m,mi,saveInverseCTF);
%     WriteJpeg(mf,[jpegPath basename '.jpg']);
%     imags(mf);
%     title(d(i).name,'interpreter','none');
%     drawnow;
% 
% %         ms=uint8(imscale(m));
% %          image(ms);
% %         imwrite(ms,[jpegPath basename '.jpg'],'jpg');
%     end;
% end;
% disp('Done.');
