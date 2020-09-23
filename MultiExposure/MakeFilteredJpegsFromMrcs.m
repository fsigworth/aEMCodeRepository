% MakeFilteredJpegsFromMrcs
% Read mi files from the Info directory. Based on that, load a merged image
% (preferably *s.mrc), inverse-CTF filter it, and save it as a jpeg in the 
% jpegPath directory.
figure(1);
clf;
disp(['Converting mrc images to jpgs in base directory ' pwd]);
allNames=f2FindInfoFiles;
targetSize=960;
suffix={'s' 'vs'};
jpegPath='jpeg_s/';
clipThreshold=5;  % include SDs.
doWrite=1;
invCTF=0.5;
invCTF=0;
finalLowpass=0.2;
CheckAndMakeDir(jpegPath,1);

for i=1:numel(allNames)
    mi=ReadMiFile(allNames{i});
    for j=1:numel(suffix) % read each row of the suffix array
        [m, mergeFullName,ok]=meReadMergedImage(mi,0,suffix{j});
        imSize=size(m);
        n=size(m,1);
        ds0=n/targetSize;
        ds=round(ds0);
        if abs(ds0-ds)>.2 % allowed tolerance in downsampling ratio
            disp(['Downsampling ratio discrepancy: ' num2str([ds0 ds]) '   ' num2str([n targetSize])]);
                return
        end;
        md=Downsample(m,round(size(m)/ds));
        mdf=GaussFilt(md,finalLowpass);
        if invCTF
           mdf=meCTFInverseFilter(mdf,mi,invCTF);
        end;
        [pa,baseName,ex]=fileparts(mergeFullName);
        jpegName=[jpegPath baseName '.jpg'];
        msc=WriteJpeg(mdf,jpegName,clipThreshold,doWrite);
        imaga(msc);
        title(jpegName,'interpreter','none');
        drawnow;
%         imwrite(ms,[jpegPath basename '.jpg'],'jpg');
    end;
end;
disp('Done.');
