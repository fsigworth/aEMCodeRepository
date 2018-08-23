% WriteTiffMultipage.m

mf=uint8(128+10*randn(4096,4096,10));
nim=size(mf,3);

t=Tiff('img.tif','w');

for img=1:nim
    cmf=mf(:,:,img);
    disp(['Writing ' num2str(img)]);
    
    ts=struct;
    ts.ImageLength=size(cmf,1);
    ts.ImageWidth=size(cmf,2);
    ts.Orientation=Tiff.Orientation.LeftBottom;
    ts.Photometric=Tiff.Photometric.Palette;
    
%     Find bounds and set the graymap
    x=cmf(:);
    n=numel(x);
    nfrac=round(.001*n);
    xs=sort(x);
    nfrac=min(n,nfrac);  % don't allow it to go beyond n
    p1=xs(nfrac);
    upperNFrac=round(n-nfrac)+1;
    upperNFrac=min(n,upperNFrac);
    p2=xs(upperNFrac);
    map=zeros(256,1);
    map(p1:p2,1)=(0:1/double(p2-p1):1)';
    map(p2:end)=1;
    ts.ColorMap=repmat(map,1,3);

    %     Set the other tags
    ts.BitsPerSample=8;
    ts.SamplesPerPixel=1;
    ts.RowsPerStrip=16;
    ts.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
    ts.Compression=Tiff.Compression.Deflate;
    ts.Software='MatlabRSC';
    ts.ImageDescription='s.compression=1; s.pixA=1.25';
    t.setTag(ts);
    t.write(cmf)
    t.writeDirectory; % create a new directory structure
end;
t.close;