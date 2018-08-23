function nTrunc=WriteZTiff(m,pixA,filename,mpars)
% function nTrunc=WriteZTiff(m,pixA,filename,mpars)
% Write a lossy compressed-tiff file from the image m.  We assume that the
% whitened, compressed image will fit into an 8-bit range; the number of
% points outside that range is returned as nTrunc.
% Version 1.01 adds the field s.size to allow getting the size from the
% tags only.
% Version 1.02 adds the user-specified field s.origName

% here are the default parameter values:
pars.origName='';  % default original filename is empty string.
pars.snrRatio=100;  % target noise-to-noise ratio
pars.lfCutoff=.07; % radius of frequency domain not fitted

pars.displayOn=0;  % disply fit and error spectra
pars.origClass=class(m);  % allows recovery of the original variable class.
pars.useDeflate=0; % 2x slower 'deflate' compression is 5% smaller

% Insert the user-specified mpars values
if nargin<4
    mpars=struct;
end;
pars=SetOptionValues(pars,mpars);

[nx,ny,nim]=size(m);

binFactor=8;  % binning for power spectrum fitting.
sFactor=pars.snrRatio/12;
n=floor([nx ny]/binFactor);  % same behavior as BinImage
t=Tiff(filename,'w');  % create the output file object
nTrunc=0;

for iimg=1:nim  % loop over the stack of images
    m1=double(m(:,:,iimg));
    fm1=fftn(m1);
    sp2=fftshift(abs(fm1).^2)/numel(m1);  % normalized spectral density
    sp2b=BinImage(sp2,binFactor);
    
    %%  Construct the prewhitening filter
    
    fs0=[.05 .07 .1 .14 .2 .4];  % Gaussian frequencies, a row vector
    fs=fs0;
    fs(fs0<pars.lfCutoff)=[];  % eliminate Gaussians that are too small
    if numel(fs)<1
        fs=fs0(end);
    end;
    mask=1-fuzzymask(n,2,n(1)*pars.lfCutoff); % no fitting below cutoff
    % We do 2 least-squares fits.  On the second round, we upweight low-valued
    % points through the use of the first fit.  This is a poor-man's way to do
    % a nonlinear fit of the log of the spectrum.
    
    weight=1;  % initial weights for LS fitting
    for j=1:2 % We fit twice
        if j>1
            weight=model.^-.5; % 2nd time is a weighted fit
        end;
        sp2m=sp2b.*mask.*weight;
        Fcns=zeros(prod(n),numel(fs));
        for i=1:numel(fs);
            f=ccdAliasedGaussians(n,fs(i),1).*mask.*weight;
            Fcns(:,i)=f(:);
        end;
        a=LinLeastSquares(Fcns,double(sp2m(:)));
        %         a=max(a,0);  % don't allow negative coefficients
        model=reshape(ccdAliasedGaussians(n,fs,a),n(1),n(2));
        model=max(model,1e-6);
    end;
    
    %     if pars.fancyDisplay
    %         ztFancyDisplay(m1,filename,sp2b,model,mask);
    %     else
    if pars.displayOn
        subplot(221);  % Show the original image
        imags(BinImage(m1,binFactor));
        [pa,nm,ex]=fileparts(filename);
        title([nm ex ' frame ' num2str(iimg)],'interpreter','none');
        axis off;
        d=(log(sp2b)-log(model)); % look at the log discrepancy
        subplot(223); imacs(d.*mask);  % show the log 2D spectrum-model
        title('2D spectral-fit error');
        % and 1D spectrum, model and masked model.
        freqs=0:1/n(1):(n(1)-1)/(2*n(1));
        subplot(224);
        plot(freqs,[Radial(log(sp2b)) Radial(log(model)) Radial(log(model).*mask)]/log(10));
        ylabel('log spectral density');
        xlabel('frequency')
        title('1D spectral fit');
        drawnow;
    end;
    
    as=a/sFactor;  % scale for output
    as=as(:)';  % Force to be a row vector for string conversion.
    fsModel=ccdAliasedGaussians(nx,fs,as);  % assume square image here!!
    filter=1./sqrt(fsModel);
    
    m2=real(ifftn(fm1.*ifftshift(filter)));
    nTrunc=sum(abs(m2(:))>127)+nTrunc;
    cmf=uint8(128+m2);  % offset binary
    
    %       Decompress here just to check fidelity
    if pars.displayOn
        rf=single(cmf)-128;
        rm=real(ifftn(fftn(rf).*ifftshift(sqrt(fsModel))));
        subplot(222);
        semilogy(freqs,[RadialPowerSpectrum(m1,1,binFactor)...
            RadialPowerSpectrum(m1-rm,1,binFactor)]);
        ylabel('log spectral density');
        xlabel('frequency')
        title([num2str(nTrunc) ' points truncated']);
        legend('Original','Compression error');
        drawnow;
    end;
    
    ts=struct;  % structure to hold the Tiff tag values
    ts.ImageLength=size(cmf,1);
    ts.ImageWidth=size(cmf,2);
    ts.Orientation=Tiff.Orientation.LeftBottom;
    ts.Photometric=Tiff.Photometric.Palette;
    
    %       Make the grayscale colormap.  We expand the few gray values to fill
    %       the 256 range, cutting 1/1000 of the histogram from each end.
    [p1,p2]=Percentile(cmf,.001);
    p1=double(max(1,min(p1,255)));
    p2=double(max(1,min(max(p1+1,p2),256)));  % Force >= 2 gray levels
    map=zeros(256,1);
    map(p1:p2,1)=(0:1/(p2-p1):1)'; % make a ramp from p1 to p2
    map(p2:end,1)=1;
    ts.ColorMap=repmat(map,1,3);  % R, G and B all the same.
    
    ts.BitsPerSample=8;
    ts.SamplesPerPixel=1;
    ts.RowsPerStrip=2^(ceil(18-log2(nx)));  % strip has max 256k pixels
    ts.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
    if pars.useDeflate
        ts.Compression=Tiff.Compression.Deflate; % Compression is set here.
    else
        ts.Compression=Tiff.Compression.LZW; % LZW is twice as fast.
    end
    ts.Software='Matlab ZTiff';
    
    %   Write the multi-line ImageDescription, which contains variable
    %   assignments in Matlab syntax.  Note: a literal string containing one
    %   quote character is ''''
    ts.ImageDescription=sprintf('%s%s%s%s;\n',...
        's.compVersion=', '',  '1.02',                '',...
        's.origName=',  '''',  pars.origName,       '''',...
        's.origPixA=',    '',  num2str(pixA),         '',...
        's.origSize=',   '[',  num2str([nx ny nim]), ']',...
        's.origClass=', '''',  pars.origClass,      '''',...
        's.size=',       '[',  num2str([nx ny nim]), ']',...
        's.snrRatio=',      '', num2str(pars.snrRatio),'',...
        's.fs=',         '[',  num2str(fs,8),        ']',...
        's.as=',         '[',  num2str(as,8),        ']');
    t.setTag(ts);
    t.write(cmf)
    if iimg<nim
        t.writeDirectory;  % In case we have more to write
    end;
end; % for iimg
t.close;
