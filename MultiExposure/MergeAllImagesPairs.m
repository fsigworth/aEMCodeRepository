% MergeAllImagesPairs.m
% special version that operates only on en and em images.

% serverDir='/Volumes/TetraData/EMWork/Liguo/FavoriteBK/';  % main directory where images are stored
% serverDir='/Volumes/TetraData/EMWork/Liguo/FavoriteBK/04T-BKPo-10dec18a-leginon-part/rawdata/';  % main directory where images are stored
serverDir='/Volumes/TetraData/EMWork/Liguo/BKfavorite2010/03T-BKPc-10oct02a-Leginon-part/0BKc-good-1002/';
imagePath='';         % subdirectory for these images
infoPath='';            % complete path to where the info files are stored.
procPath=imagePath;
mcDS=2;             % Composite image downsample ratio
kWiener=.04;        % Wiener filter constant
fc=.1;              % Lowpass filter corner frequency

% load or, if necessary, create the micrograph info structures
overwrite=0
infos=meScanLeginonFiles2(serverDir,imagePath,procPath,infoPath,overwrite);
nmi=numel(infos);

%%
for i=1:nmi
    
    mi=infos{i};  % pick up the structure
    mi.imagePath=imagePath;
    mi.procPath=procPath;
    mi.infopath=infoPath;
    
    % Construct the complete image file names
    nim=numel(mi.rawFile);
    if nim>2
        mi.rawFile=mi.rawFile(1:2);
        nim=2;
    end;
    fnames=cell(1,nim);
    disp('Processing images:');
    for i=1:nim
        fnames{i}=[serverDir mi.imagePath char(mi.rawFile{i})];
        disp(fnames{i});
    end;
    
    % Read the images
    [m mi.pixA mi.doses]=meReadImages(fnames);
    sz=size(m);
    mi.nPix=sz(1:2);
    
    % operate with the pre-whitening filter
    disp('mePreWhiten');
    m=mePreWhiten(m);
    
    % Fit the ctfs
    mi.ctf=meFitCTFs(m,mi.pixA,[1.5 4 15],1);
 %%   
    % Align the images
    mi.mergeMatrix=meAlignMultiExposures(m,mi.pixA,mi.ctf,mi.doses,4);
%%    
    % Store the CTF and alignment information
%     save([serverDir infoPath mi.baseFilename 'mi2.mat'],'mi');

    % Create the composite image mc
    [effctf mc]=meCombineImages(m,mi.pixA,mi.ctf,mi.mergeMatrix,mi.doses,mcDS);
    mc=-mc;  % go back to reversed contrast
    
%     [effctf mc]=rsCombineImages(m,mi,mcDS);
    %%
    % Filter the combined image
    mcf=meFilterImage(mc,effctf,mcDS*mi.pixA,kWiener,fc);  % Filter the combined image
    
    % Write out the combined image and filtered version
    basename=[serverDir mi.imagePath mi.baseFilename];
    WriteMRC(mc,mcDS*mi.pixA,[basename 'm2.mrc']);
    WriteMRC(mcf,mcDS*mi.pixA,[basename 'mf2.mrc']);
    imwrite(uint8(imscale(mcf)),[basename 'mf2.jpg']);
    imacs(mcf);
    title(basename,'interpreter','none');
%%
    mcs=meRemoveCrystal(mc,mcDS*mi.pixA,6,1);  % crystal-subtracted image
    WriteMRC(mcs,mcDS*mi.pixA,[basename 'mc2.mrc']);
    imwrite(uint8(imscale(mcs)),[basename 'mc2.tif']);

 %%   
    
 
 mi=meFindVesicles(mcs, mi);
 
 disp('Reconstructing vesicles');
 vm=meMakeModelVesicles(mc, mi);
 mcv=mcs-vm;
 %%
 subplot(2,3,1);
 imacs(BinImage(mcs,4));
 subplot(2,3,4);
 imacs(BinImage(mcv,4));
 title('Subtracted');
 drawnow;
 
    save([serverDir infoPath mi.baseFilename 'mi2.mat'],'mi');
    WriteMRC(mcv,mcDS*mi.pixA,[basename 'mcv2.mrc']);
    imwrite(uint8(imscale(mcv)),[basename 'mcv2.tif']);
disp(['TIFF pixel size is ' num2str(mi.pixA*mcDS) ' A']);
    



end;
toc
nmi
