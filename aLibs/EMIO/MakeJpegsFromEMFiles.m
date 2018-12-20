% function MakeJpegsFromEMFiles(OutputDir, binning, display, defaultPixA)
% function MakeJpegsFromEMFiles(binning, display)
% For each EM image file (.mrc, .dm3, .img, .hed, .tif)
% in the current directory, make a .jpg file.  All arguments are optional.
% If OutputDir is given and is not an empty string, the
% jpeg files are written there.  The argument binning controls the degree
% of downsampling.  Its default value is 2.  If display=0 then no display
% is shown. DefaultPixA controls the labeling of the display.
% fs July 2009 rev. Apr 2011, Nov 2017

% inputExtensions={'.mrc' '.tif' '.dm4' '.mrcs'};
inputExtensions={'.mrc' '.tif' '.mrcs'};
% inputExtensions={'.mrcs'};
sumStacks=1;

outputDir='/ysm-gpfs/project/fjs2/181216/No5Graphene/sq05_1/JpegSums/';

disp(['Converting EM files to jpgs from directory ' pwd]);
nargin=0;

% if nargin<1
%     OutputDir='';
% end;
len=numel(outputDir);
if len>0
    if outputDir(len)~='/' % doesn't end with slash
        outputDir=[outputDir '/'];
    end;
    disp(['Writing output files to ' outputDir]);
end;
if nargin<2
    binning=4;
end;
if nargin<3
    display=1;
end;
if nargin<4
    defaultPixA=0;
    defaultPixA=1.87;
end;
d=dir;
for i=3:numel(d)
    disp(d(i).name);
    [~,~,ex]=fileparts(d(i).name);
    ok=any(strcmp(ex,inputExtensions)); % one of the image file types
    if ok
        if strcmp(ex,'.mrcs')
            [mi, s, ok]=ReadMovie(d(i).name);
            pixA=s.pixA;
        else
            [mi, pixA, ok]=ReadEMFile(d(i).name);
        end;
    else
        continue;
    end;
    nim=size(mi,3)
    ok=ok && (sumStacks ||size(mi,3)==1);  % don't do stacks
    if ok
        disp(d(i).name);
        if nim>1
            disp(['Summing ' num2str(nim) ' frames.']);
            mi=sum(single(mi),3);
        end;
        m=RemoveOutliers(mi,4);
        n=size(m,1);
        n=NextNiceNumber(n,5,4);
        me=mean(m(:));
        m=Crop(m,n,1,me);
        if binning>1
            m=Downsample(m,n/binning);
        end;
        ms=uint8(imscale(m,256,1e-4));
        if display
            n=size(ms);
            pixA=max(pixA,defaultPixA)*binning;
            if pixA==0
                pixA=1;
                label='pixels';
            elseif pixA*n(1)>1000
                pixA=pixA/10;
                label='nm';
            else
                label='A';
            end;
            figure(1); SetGrayscale;
            imaga((1:n(1))*pixA,(1:n(2))*pixA,ms);
            title([d(i).name '    original pixel size = ' num2str(pixA/binning) label],'interpreter','none');
            xlabel(['Dimension, ' label]);
            drawnow;
        end;
        [~, nm]=fileparts(d(i).name);
        imwrite(ms,[outputDir nm '.jpg'],'jpg');
    else
        disp(['Skipped: ' d(i).name]);
    end;
end;
disp('Done.');
