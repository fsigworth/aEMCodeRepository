function f2Compressor(arg1,arg2)

binFactor=8;
snrRatio=50;
% snrRatio=MyInput('snrRatio',50);

figure(1);
set(gcf,'name','f2Compressor');

[names, pathName]=uigetfile('*.*','Select image files','multiselect','on');
if isnumeric(pathName) % File selection was cancelled
    return
end;
if ~iscell(names)
    names={names};
end;
cd(pathName);


nFiles=numel(names);
for fIndex=1:nFiles
    name=names{fIndex};
    disp(['Reading ' name]);
    if ~strcmp(name(end-5:end),'.z.tif')  % not already compressed
        [m,pixA,ok]=ReadEMFile(name);
    else
        ok=0;
    end;
    if ok % Valid file
%         tic
        mf=single(m);

        [nx,ny,nim]=size(mf);
%         toc
        
        [pa,nm,ext]=fileparts(name);
        outName=[nm '.z.tif'];
        disp(['Output name: ' outName]);
        t=Tiff(outName,'w');  % create the output name
        
        sFactor=snrRatio/12;
        n=[nx ny]/binFactor;
%         disp('Outliers');
%         tic
        md1=RemoveOutliers(mf);
%         toc;
%         disp('Compress');
%         tic
%         imags(md1)
        fmd1=fftn(md1);
        sp2=fftshift(abs(fmd1).^2)/numel(md1);
        sp2b=BinImage(sp2,binFactor);
        
        %%  Construct the prewhitening filter
        
        % na=1;
        fs=[.07 .1 .14 .2 .4];
        %     mask=1;
        mask=1-fuzzymask(n,2,n(1)/20,n(1)/50);
        weight=1;
        for j=1:2
            if j>1
                weight=model.^-.5;
            end;
            % weight=1./sqrt(model);
            
            sp2m=sp2b.*mask.*weight;
            
            Fcns=zeros(prod(n),numel(fs));
            for i=1:numel(fs);
                f=ccdAliasedGaussians(n,fs(i),1).*mask.*weight;
                %     f=amSpectrum2(n,fs(i),1,na).*mask.*weight;
                Fcns(:,i)=f(:);
            end;
            a=LinLeastSquares(Fcns,double(sp2m(:)));
            model=reshape(ccdAliasedGaussians(n,fs,a),n(1),n(2));
            model(model<1e-6)=1e-6;
        end;
        
        d=(log(sp2b)-log(model));
%         subplot(221); imacs(d);
%         subplot(222); plot(d); %axis([-inf inf -1 .5])
        clf;
        subplot(221);
        imags(BinImage(md1,binFactor));
        title(outName,'interpreter','none');
        axis off;
        subplot(223); imacs(d.*mask);
        subplot(224); plot([Radial(log(sp2b)) Radial(log(model)) Radial(log(model).*mask)]/log(10));
        drawnow;
        
        as=a/sFactor;  % scale for output
        fsModel=ccdAliasedGaussians(nx,fs,as);  % assume square image
        filter=1./sqrt(fsModel);
        
%         plot(sect(filter));
%         drawnow;
        %
        mf=real(ifftn(fmd1.*ifftshift(filter)));
        cmf=uint8(128+mf);
        % cmf=uint16(128+mf);
%         toc;
         disp('Writing');
%         tic

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
%         ts.RowsPerStrip=16;
        ts.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
        ts.Compression=Tiff.Compression.Deflate;
        ts.Software='MatlabRSC';
        parString=['s.fs=[' sprintf('%g ',fs) ']; s.as=[' sprintf('%g ',as) '];'];
        ts.ImageDescription=['s.compression=1; s.pixA=' num2str(pixA) '; ' parString];
        %     ts.TransferFunction=(1:256);  % apparently uint16s.
        % ts.ModelTransformationMatrixTag=zeros(1,16);
        t.setTag(ts);
        % vals=zeros(1,10,'uint16');
        % t.setTag(65011,vals);
        t.write(cmf)
        t.close;
%         toc
        %         Show compression performance
%         disp('Plotting spectra');
%         tic
        bf=16;
        rf=single(cmf-128);
        rm=real(ifftn(fftn(rf).*ifftshift(sqrt(fsModel))));
        subplot(222);
        semilogy([RadialPowerSpectrum(md1,1,bf) RadialPowerSpectrum(md1-rm,1,bf)]);
        drawnow;
        disp(' done.');
%         toc
    end; % if ok
end; % for fIndex

pause(2);
close(1);
