% meWriteJpegs
% From a set of images, create
% jpeg copies of them.  Store them in a ./Jpeg folder.

jpegPath='Jpeg/';
gaussFc=.1;
outputBinning=2;
clipThreshold=1e-3;

[fname, imagePath]=uigetfile('*.*','Select image files','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;
if isnumeric(imagePath)  % user clicked Cancel
    return
end;
cd(imagePath);
%%
figure(1);
SetGrayscale;

nfi=numel(fname);
for fileIndex=1:nfi
    name=fname{fileIndex};
    [m, pixA, ok]=ReadEMFile(name);
    if ok
        m=single(m);
        if gaussFc>0
            m=GaussFilt(m,gaussFc);
        end;
        if outputBinning>1
            m=BinImage(m,outputBinning);
        end;
        m=imscale(m,256,clipThreshold);
        imacs(m);
        title(name,'interpreter','none');
        drawnow;
        if ~exist(jpegPath,'dir')
            mkdir(jpegPath);
        end;
        [pa, baseName]=fileparts(name);
        jpegName=[jpegPath baseName '.jpg'];
        disp(jpegName);
        imwrite(uint8(rot90(m)),jpegName);  % WriteJpeg(m,jpegName);
    else
        disp(name);  % warning message is already displayed by ReadEMFile
    end;
end;
