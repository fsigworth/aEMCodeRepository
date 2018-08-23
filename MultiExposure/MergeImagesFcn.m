function [mi,ok,effCTF] = MergeImagesFcn(mi,pars,ctfOptions,modelSpectrum)
% the structure pars has fields
% overwrite
% doFitting
% doAlignment
% ignoreOldCTFs
% initialDefoci
% ctfJPath
% mcDS
% nZeros
% useCircMask
% mergeMode
% suffix
% logs
% useGraphics
% writeGraphics
% overridePixA
%
ctfOptions.showGraphics=pars.showGraphics;  % we set only this ctfOption, and kV.
ctfOptions.kV=mi.kV;

% Check to make sure we have things to do
mergeDataPresent=numel(mi.mergeMatrix)>3;
ok=false;
if numel(mi.imageFilenames)<1
    mdisp(pars.logs,['No image files: ' mi.baseFilename]);
    return
elseif mergeDataPresent && ~pars.overwrite
    mdisp(pars.logs,['Already merged: ' mi.baseFilename]);% alignment has already been done
    return
else
    doFit=(~mergeDataPresent || pars.doFitting);
    doAlign=(~mergeDataPresent || pars.doAlignment);
end;
ok=true;

% Read the images
if ischar(mi.imageFilenames) % Make sure this is a cell array.
    mi.imageFilenames={mi.imageFilenames};
end;
nim=numel(mi.imageFilenames);
if numel(mi.weights)<nim
    mi.weights(end:nim)=1;
elseif numel(mi.weights)>nim
    mi.weights=mi.weights(1:nim);
end;
[m, tempPixA, doses]=meReadImagesNorm(mi,mi.cpe,0,pars.overridePixA,mi.weights,1);
if ~isfield(mi,'doses') || (~any(isnan(doses)) && any(doses>.1))
    mi.doses=doses;
end;
if tempPixA>0.1 && (tempPixA~=mi.pixA)
    mi.pixA=tempPixA;
    disp(['Setting mi.pixA to ' num2str(mi.pixA)]);    
end;

%nim=min(size(m,3),size(mi.frameSets,1));
%m=m(:,:,1:nim);
nim=size(m,3);
mi.doses=mi.doses(1:nim);


if pars.removeOutliers
    for i=1:nim
        m(:,:,i)=RemoveOutliers(m(:,:,i));
    end;
end;

if numel(pars.weights)<nim
    pars.weights=single(ones(1,nim));
end;
mi.weights=pars.weights;

mdisp(pars.logs,['pixA = ' num2str(mi.pixA)]);
mdisp(pars.logs,['doses = ' num2str(mi.doses.*mi.weights)]);
sz=size(m);
mi.imageSize=sz(1:2);
modelSpectrum=Crop(modelSpectrum,mi.imageSize);  % possible downsampling
% operate with the pre-whitening filter
%     mdisp(pars.logs,'mePreWhiten');
m=mePreWhiten(m,modelSpectrum);
%%
if doFit
    % Fit the ctfs
    %         ncts=numel(mi.ctf);
    rawCTFs=struct;
    defoci=zeros(2,nim);
    for j=nim:-1:1 % ------actual fitting is done here------
        [rawCTF,~,ctfDisDat]=meFitCTF(m(:,:,j),mi,pars,ctfOptions,j,0); 
        ctfDisDat.contour.title=mi.imageFilenames{j};
        if j==nim
            rawCTFs=rawCTF;
        end;
        rawCTFs(j)=rawCTF;
        defoci(:,j)=[rawCTFs(j).defocus; rawCTFs(j).deltadef];
        if pars.showGraphics && pars.writeGraphics
            jctfName=[mi.basePath pars.ctfJPath mi.baseFilename num2str(j) '-spect.jpg'];
            print('-djpeg','-r150',jctfName);  % save the CTF window.
        else
           ctfDatName=[mi.basePath pars.ctfJPath mi.baseFilename 'al' char(96+j) '-ctfDisDat.mat'];
           save(ctfDatName,'ctfDisDat');
        end;
    end;
    mdisp(pars.logs,['defocus = ' num2str(defoci(1,:))]);
    mdisp(pars.logs,['astig   = ' num2str(defoci(2,:))]);
%     subplot(2,3,2);  % use this port for subsequent image display
    %%
    % Align the images
    % use downsampling by 4
    if doAlign  % we'll skip the alignment step.
        oldMats=[];
    else
        oldMats=mi.mergeMatrix;
    end;
        
    [mi.mergeMatrix, mi.ctf, mi.mergeSNR]=meAlignMultiExposures(m,mi.pixA,rawCTFs,mi.doses,4,oldMats,pars.showGraphics);
    %             Ps
    %%
    % Store the CTF and alignment information
else
    mdisp(pars.logs,'Fitting skipped.');
end;
%%
% Create the composite image mc
[effCTF, mc]=meCombineImages2(m,mi,pars.mcDS,pars.useCircMask,pars.nZeros,pars.mergeMode);
mi.mergeVersion=3;
mi.mergeMode=pars.mergeMode;
% Note that the effctf doesn't include ccd transfer function, but
% we're not using it anyway for merging.
if pars.mergeMode~=3
    mc=-mc;  % go back to reversed contrast
end;
%             Make a mask indicating the valid edge of the image.
mskSize=round(mi.imageSize(1)/4);
msk=meMakeMergedImageMask(mskSize,mi.mergeMatrix);
mi=meInsertMask(msk,mi,1);  % merge mask is the first one in the mask stack.
%%
% Write out the combined image and filtered version.
procname=[mi.basePath mi.procPath mi.baseFilename];  % Use the processed image directory
jpegname=[mi.basePath pars.mergeJPath mi.baseFilename];

% For jpeg files we trim .1% of values from the intensity
% histograms.
if pars.writeZTiff
    ztPars.snrRatio=200;  % target noise-to-noise ratio
    ztPars.lfCutoff=.2; % radius of frequency domain not fitted
    mergeOutName=[procname pars.suffix 'z.tif'];
    WriteZTiff(mc,pars.mcDS*mi.pixA,mergeOutName,ztPars);
    mdisp(pars.logs,['Wrote merged image: ' mergeOutName]);
end;
if pars.writeZTiff<2
    mergeOutName=[procname pars.suffix '.mrc'];
    WriteMRC(mc,pars.mcDS*mi.pixA,mergeOutName);
    mdisp(pars.logs,['Wrote merged image: ' mergeOutName]);
    WriteMRC(Downsample(mc,pars.smallImageSize),mi.pixA*mi.imageSize(1)/pars.smallImageSize, ...
            [procname pars.suffix 's.mrc']);
end;
WriteJpeg(BinImage(mc,2),[jpegname pars.suffix '.jpg']);

% Show the image and spectrum
if pars.showGraphics
     subplot(2,3,2);
    imacs(BinImage(mc,4));
    axis off;
    title(mergeOutName,'interpreter','none');
    if pars.writeGraphics
%         figure(10)
        QuickLookSpectrum(mc,pars.mcDS*mi.pixA,[mi.baseFilename pars.suffix '.mrc']);
        jSpectName=[mi.basePath pars.ctfJPath mi.baseFilename 'cs.jpg'];
        print('-djpeg','-r150',jSpectName);  % save the CTF window.
    end;
end;