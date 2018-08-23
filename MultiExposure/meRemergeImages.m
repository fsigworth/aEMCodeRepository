% meRemergeImages.m
%  **maybe not valid for version > 12 **
% Given mi files, compute the merged images using weights to specify which
% images to use.
% We save the mi, m, jpg files in new directories, e.g. Info_110 and Merged_110
% fs 16 Apr 13

cpe=16;
mcDS=2;  % Downsampling of merged image
iCamera=1;  % F20 ultrascan
weights=[1 1 1];
makeDirectories=1;  % Create new directories if needed

weightString=sprintf('%1d',weights);
modelSpectrum=CCDModelSpectrum2D(iCamera);  % handle DE-12 or CCD

% Ask for a set of original mi files.
[fname pa]=uigetfile('*mi.mat','Select mi files to duplicate','multiselect','on');
[rootPath infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;

cd(rootPath);
nfiles=numel(fname);
figure(1);
clf;
SetGrayscale;

%%
for findex=1:nfiles
    % Get the mi structure
    load([infoPath fname{findex}]);
    if ~isfield(mi,'tempPath')
        mi.tempPath='';
    end;
    
    mi0=mi;  % keep the old one
%     update the weights
    mi.weights=weights;
    
    % Create directories
    jpegPath=[mi0.procPath(1:end-1) '/jpeg/'];
% 
%     mi.infoPath=[mi0.infoPath(1:numel(mi0.infoPath)-1) '_' weightString '/'];
%     mi.procPath=[mi0.procPath(1:numel(mi0.procPath)-1) '_' weightString '/'];
%     jpegPath=[mi0.procPath(1:numel(mi0.procPath)-1) '_' weightString '/jpeg/'];
%     mi.tempPath=['Temp_' weightString '/'];
    if findex==1
        writePaths={mi.infoPath; mi.procPath; jpegPath; mi.tempPath};
        %if needed, create directories.  The imagePath must already exist.
        for j=1:numel(writePaths);
            if ~exist(writePaths{j},'dir')
                if makeDirectories
                    mkdir(writePaths{j});
                    disp([writePaths{j} ' created']);
                else
                    error(['the path doesn''t exist: ' writePaths{j}]);
                end;
            end;
        end;
    end;

    % Read the images
    nim=numel(mi.imageFilenames);
    fnames=cell(1,nim);
    disp('Processing images:');
    for j=1:nim
        fnames{j}=[mi.imagePath mi.imageFilenames{j}];
    end;
    
    m=meReadImages(fnames,cpe,0,mi.pixA,weights);
    m=mePreWhiten(m,modelSpectrum);
    
    %%
    if weights(1)==0
        ds=1;
    else
        ds=mcDS;
    end;
    [effctf mc mts]=meCombineImages(m,mi,ds);
    mc=-Downsample(mc,mi.imageSize/mcDS);   % go back to reversed contrast
    %             Make a mask indicating the valid edge of the image.
    mskSize=round(mi.imageSize(1)/4);
    msk=meMakeMergedImageMask(mskSize,mi.mergeMatrix, 0, weights);
    mi=meInsertMask(msk,mi,1);  % merge mask is the first one in the mask stack.
    % Write out the combined image.
    % For jpeg files we trim .1% of values from the intensity
    % histograms.
    WriteMRC(mc,mcDS*mi.pixA,[mi.procPath mi.baseFilename 'm.mrc']);
    WriteJpeg(mc,[jpegPath mi.baseFilename 'm.jpg']);
    disp(['Wrote merged images: ' mi.procPath mi.baseFilename 'm.mrc']);
    figure(1);
    imacs(mc);
    title(mi.baseFilename,'interpreter','none');
    
    % Show the image and spectrum
    figure(10);
    QuickLookSpectrum(mc,mcDS*mi.pixA,[mi.baseFilename 'm.mrc']);
    save([mi.infoPath fname{findex}],'mi');
    
end;
