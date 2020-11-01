% MakeFilteredJpegsFromMrcs
% Read mi files from the Info directory. Based on that, load a merged image
% (preferably *s.mrc), inverse-CTF filter it, and save it as a jpeg in the
% jpegPath directory.
figure(1);
clf;
disp(['Converting mrc images to jpgs in base directory ' pwd]);

infoDir='Info';

allNames=f2FindInfoFiles(infoDir);
targetSize=768;
targetSize=1024;
% suffix={'s' 'vs'};
suffix={'s'};
readZTiff=0;
jpegPath='jpeg_s/';
clipThreshold=5;  % include SDs.

doWrite=1;
invCTF=0.3;
finalLowpass=0.2;

maxImages=200;
iStart=1;

CheckAndMakeDir(jpegPath,1);

ni=min(numel(allNames)-iStart+1,maxImages);
disp(['Working on ' num2str(ni) ' images.']);


for i=iStart:iStart-1+ni
    mi=ReadMiFile(allNames{i});
    for j=1:numel(suffix) % read each row of the suffix array
        if readZTiff
            [m,pixA,ok]=ReadEMFile(['Merged/' mi.baseFilename 'mz.tif']);
        else
            [m, mergeFullName,ok]=meReadMergedImage(mi,0,suffix{j});
        end;
        if ~ok
            continue;
        end;
        imSize=size(m);
        n=size(m,1);
        ds0=n/targetSize;
        ds=round(ds0);
        if abs(1-ds/ds0)>.2 % allowed tolerance in downsampling ratio
            disp(['Downsampling ratio discrepancy: ' num2str([ds0 ds]) '   ' num2str([n targetSize])]);
            return
        end;
        md=Downsample(m,round(size(m)/ds));
        mdf=GaussFilt(md,finalLowpass);
        if invCTF
            mdf=meCTFInverseFilter(mdf,mi,invCTF);
        end;
        defStr=[num2str(mi.ctf(1).defocus,3) 'um'];
        invStr=[num2str(invCTF,2) 'inv'];
        extraString=['  ' defStr ' ' invStr];
        if readZTiff
            baseName=mi.baseFilename;
        else
            [pa,baseName,ex]=fileparts(mergeFullName);
        end;
        
        jpegName=[jpegPath baseName '_' defStr '_' invStr '.jpg'];
        msc=WriteJpeg(mdf,jpegName,clipThreshold,doWrite);
        imaga(msc);
        title([num2str(i) '  ' jpegName],'interpreter','none');
        drawnow;
    end;
end;
disp('Done.');
