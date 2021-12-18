function MakeJpegsFromEMFiles(pars)
% function MakeJpegsFromEMFiles(pars)
% For each EM image file (.mrc, .tif, .mrcs)
% in the current directory, make a graphics file, by default a .jpg file,
% but pars.outputExtension can specify .tif or .png instead.
% If OutputDir is given and is not an empty string, the
%  files are written there.  The argument binning controls the degree
% of downsampling.  Its default value is 1.  If doDisplay=0 then no display
% is shown. DefaultPixA affects only the labeling of the display.
% fs July 2009 rev. Apr 2011, Nov 2017, Nov 2021.
if nargin<1
    pars=struct;
end;

dpars.inputExtensions={'.mrc' '.tif' '.mrcs'};
dpars.suffixes={'ms' 'mvs' '_u' '_v'}; % {} ignores filename suffix

dpars.outputExtension='.jpg';
dpars.sumStacks=1; % sum the frames in an .mrcs file
dpars.outputDir='Jpeg/'; % inside the current dir
dpars.binning=4;
dpars.doDisplay=1;
dpars.defaultPixA=0; % used for local display
dpars.clipThreshold=1e-3;  % saturation of black and white pixels
dpars.overwrite=0;
dpars.finalIndex=inf;  % write every file in the current directory
dpars.listFilenames=1;  % show the input filenames

pars=SetDefaultValues(dpars,pars,1);


% outputDir='/ysm-gpfs/project/fjs2/181216/No5Graphene/sq05_1/JpegSums/';

disp(['Converting EM files to ' pars.outputExtension 's from directory ' pwd]);

CheckAndMakeDir(pars.outputDir,1);

disp(['Writing output files to ' pars.outputDir]);

d=dir;
index=0;

for i=1:numel(d);
    [~,nm,ex]=fileparts(d(i).name);
    ok=any(strcmp(ex,pars.inputExtensions)); % one of the image file types
    sufOk=numel(pars.suffixes)==0;
    for j=1:numel(pars.suffixes)
        sufOk=sufOk || strndcmp(nm,pars.suffixes{j});
    end;
    if ok && sufOk
        if pars.listFilenames
            disp(d(i).name);
        end;
        if strcmp(ex,'.mrcs')
            [mi, s, ok]=ReadMovie(d(i).name);
            pixA=s.pixA;
        else
            [mi, pixA, ok]=ReadEMFile(d(i).name);
%             Code for writing out vesicle models. We read *mvs files,
%        nm(end-1)=[];   % and convert the name
%         miv=ReadEMFile([nm ex]);  % to read the *ms file
%         mi=miv-mi;  % The difference is written out.

% %         code for writing out vesicle center lines as *mvs
% %         mif=ReadMiFile(['../Info/' nm(1:end-2) 'mi.txt']);
% %         mif.vesicleModel=zeros(49,1);
% %         mif.vesicleModel(23:27)=1;
% %         mif.vesicle.extraPeaks=[];
% %         mi=-meMakeModelVesicles(mif,960,0,0,0,1);
% %         imags(mi);

        end;
    else
        continue;
    end;
    nim=size(mi,3);
    ok=ok && (pars.sumStacks || size(mi,3)==1);  % don't do stacks
    if ok && index<pars.finalIndex
        index=index+1;
%         disp(d(i).name);
        if nim>1
            disp(['Summing ' num2str(nim) ' frames.']);
            mi=sum(single(mi),3);
        end;
        m=RemoveOutliers(mi,4);
        n=size(m,1);
        n=NextNiceNumber(n,5,4);
        me=mean(m(:));
        m=Crop(m,n,1,me);
        if pars.binning>1
            m=Downsample(m,n/pars.binning);
        end;
        ms=uint8(imscale(m,256,pars.clipThreshold));
        if pars.doDisplay
            n=size(ms);
            pixA=max(pixA,pars.defaultPixA)*pars.binning;
            if pixA==0
                pixA=1;
                label='pixels';
            elseif pixA*n(1)>1e4
                pixA=pixA/10;
                label='nm';
            else
                label='A';
            end;
            figure(1); SetGrayscale;
            imaga((1:n(1))*pixA,(1:n(2))*pixA,ms);
            title([d(i).name '    original pixel size = ' num2str(pixA/pars.binning) label],'interpreter','none');
            xlabel(['Dimension, ' label]);
            drawnow;
        end;
        [~, nm]=fileparts(d(i).name);
        fullOutName=[pars.outputDir nm pars.outputExtension];
        if ~exist(fullOutName,'file') || pars.overwrite
            %             imwrite(ms,fullOutName,'jpg');
            WriteJpeg(m,fullOutName,0);
        else
            disp(' --exists, not written.');
        end;
    elseif index >= pars.finalIndex
        break;
%         disp(['Skipped: ' d(i).name]);
    end;
end;
disp('Done.');
