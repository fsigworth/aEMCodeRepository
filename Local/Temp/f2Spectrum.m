% f2Spectrum.m
binFactor=8;

disp('Reading');
tic
cd /Users/fred/EMWork/Hideki/141220P/movie_frames/sq02
name='Dec21_17.16.44.mrc';
[m,mrcs]=ReadMovie(name);
mf=single(m);
nim=size(mf,3);
toc

[pa,nm,ext]=fileparts(name);
outName=[nm 'p.tif'];
t=Tiff(outName,'w');

nim=2;

t.WriteDirectory(nim);

snrRatio=50;
    sFactor=sqrt(snrRatio/12);
    n=4096/binFactor;
disp('Compressing and writing');
tic
for img=1:nim
    img=1;
    md=mf(:,:,img);
    disp('Outliers');
    tic
    md1=RemoveOutliers(md);
    toc;
    disp('Compress');
    tic
    imags(md1)
    fmd1=fftn(md1);
    sp2=fftshift(abs(fmd1).^2)/numel(md1);
    sp2b=BinImage(sp2,binFactor);
    
    %%  Construct the prewhitening filter
    
    figure(1);
    % na=1;
    fs=[.07 .1 .14 .2 .4];
%     mask=1;
    mask=fuzzymask(n,2,n*100,2)-fuzzymask(n,2,n/20,n/50);
    weight=1;
    for j=1:2
        if j>1
            weight=model.^-.5;
        end;
        % weight=1./sqrt(model);
        
        sp2m=sp2b.*mask.*weight;
        
        Fcns=zeros(n^2,numel(fs));
        for i=1:numel(fs);
            f=ccdAliasedGaussians(n,fs(i),1).*mask.*weight;
            %     f=amSpectrum2(n,fs(i),1,na).*mask.*weight;
            Fcns(:,i)=f(:);
        end;
        as=LinLeastSquares(Fcns,double(sp2m(:)));
        model=reshape(ccdAliasedGaussians(n,fs,as),n,n);
        model(model<1e-6)=1e-6;
    end;
    
    d=(log(sp2b)-log(model));
    subplot(221); imacs(d);
    subplot(222); plot(d); %axis([-inf inf -1 .5])
    subplot(223); imacs(d.*mask);
    subplot(224); plot([Radial(log(sp2b)) Radial(log(model)) Radial(log(model).*mask)]/log(10));
    drawnow;
    
    fsModel=ccdAliasedGaussians(4096,fs,as);
    filter=1./sqrt(fsModel);
    subplot(221);
    plot(sect(filter));
    drawnow;
    %
    mf=sFactor*real(ifftn(fmd1.*ifftshift(filter)));
    cmf=uint8(128+mf);
    % cmf=uint16(128+mf);
    toc;
    disp('Writing');
    
    ts=struct;
%     if img>nim
%         t.writeDirectory;  % create another IFD
%     else
%         ts.IFD=nim-1;
%     end;    
    ts.ImageLength=size(cmf,1);
    ts.ImageWidth=size(cmf,2);
    ts.Orientation=Tiff.Orientation.LeftBottom;
    ts.Photometric=Tiff.Photometric.Palette;
    [p1,p2]=Percentile(cmf,.001);
    if p2<=p1
        p2=p1+1;
    end;
    map=zeros(256,1);
    map(p1:p2,1)=(0:1/double(p2-p1):1)';
    map(p2:end)=1;
    ts.ColorMap=repmat(map,1,3);
    ts.BitsPerSample=8;
    ts.SamplesPerPixel=1;
    ts.RowsPerStrip=16;
    ts.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
    ts.Compression=Tiff.Compression.Deflate;
    ts.Software='MatlabRSC';
    parString=['fs=[' sprintf('%g ',fs) ']; as=[' sprintf('%g ',as) '];'];
    ts.ImageDescription=['s.compression=1; s.pixA=' num2str(mrcs.pixA) '; ' parString];
%     ts.TransferFunction=(1:256);  % apparently uint16s.
    % ts.ModelTransformationMatrixTag=zeros(1,16);
    t.setTag(ts);
    % vals=zeros(1,10,'uint16');
    % t.setTag(65011,vals);
    t.write(cmf)
end;

t.close;
% imwrite(cmf,'img.tif','Compression','deflate');
toc;
%%
disp('Reading');
tic
crf=imread('img.tif');
toc;
disp('Decompress');
tic
rf=single(crf-128);
rm=1/sFactor*real(ifftn(fftn(rf).*ifftshift(1./filter)));
rmr=round(rm);
toc
%
figure(2);
subplot(221); imags(GaussFilt(md1,.1));
subplot(222); imags(GaussFilt(rmr,.1));
subplot(224);
semilogy([RadialPowerSpectrum(md1) RadialPowerSpectrum(md1-rm) RadialPowerSpectrum(md1-rmr)]);

% subplot(223);
% plot([sect(mf) sect(round(mf)-mf)]);
subplot(223);
plot([sect(rm) sect(md1-rm)+300]);
