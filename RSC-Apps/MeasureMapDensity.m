% MeasureMapDensity.m

ds=1;  % binning factor
ok=true;

while ok
    % ----------------Get the stack files---------------
    
    [fname, pa]=uigetfile('*.mdoc','Select a map doc file');
    if isnumeric(pa)  % user clicked Cancel
        return
    end;
    
    cd(pa);
    
    %%
    stackName=fname(1:end-5);
    disp(['Reading ' stackName]);
    [st, s]=ReadMRC(stackName);
    pixA=s.pixA;
    [nx, ny, nz]=size(st);
    disp([nx ny nz]);
    
    nxw=2*floor(nx/ds/2);
    nyw=2*floor(ny/ds/2);
    disp('Downsampling');
    % the image size had better be a multiple of ds
    isStack=nz>1;
    fst=Downsample(Crop(single(st),[nxw nyw]*ds,isStack),[nxw nyw],isStack);
    
    mid=mean(fst(:));
    
    
    disp(['reading ' fname]);
    %%
    mdoc=fopen(fname);
    lines={};
    i=0;
    j=0;
    zVals=[];
    assigned=false(0,2);
    s=fgetl(mdoc);
    
    % Get image info
    mosaic=zeros(0,1,'single');
    pixelSpacing=0;
    while ~feof(mdoc) && numel(s)>0
        c=textscan(s,'%s = %f %f','emptyvalue',0);
        tok=char(c{1});
        switch tok
            case 'PixelSpacing'
                pixelSpacing=c{2};
        end;
        s=fgetl(mdoc);
    end;
    if pixelSpacing==0
        error(['No header in ' fname]);
    end;
    
    %  a blank line follows the header.
    coords=zeros(0,2);
    stage=zeros(0,2);
    while ~feof(mdoc)
        if numel(s)>0
            c=textscan(s,'%s = %f %f %f','emptyvalue',0);
            tok=char(c{1});
            switch tok
                case '[ZValue'
                    iz=c{2}+1;
                case 'PieceCoordinates'
                    xpos=floor(c{2}/ds);
                    ypos=floor(c{3}/ds);
                    coords(iz,:)=[c{2} c{3}];
                    mosaic(xpos+1:xpos+nxw,ypos+1:ypos+nyw)=fst(:,:,iz);
                    assigned(xpos+1:xpos+nxw,ypos+1:ypos+nyw)=true;
                case 'StagePosition'
                    stage(iz,:)=[c{2} c{3}];
            end;
        end;
        s=fgetl(mdoc);
    end;
    if size(coords,1)<1  % not a mosaic, just get the first image
        disp('Not a mosaic, reading one image');
        coords(1,:)=[1 1];
        mosaic=fst(:,:,1);
        assigned=true(nxw,nyw);
    end;
    % mosaic(~assigned)=refIntensity*exp(-thicknessRange*0.7/iceConstant);
    mosaic(~assigned)=0;
    
    %%
    % norm. image intensity = exp(-thk/iceConstant)
    % expand to 256 when norm intensity=0, 0 when
    % normIntensity=exp(-thickRange/iceConstant).
    % hence relthk=(-iceConstant*log(intensity))/thicknessRange;
    % mosaic=GaussFiltDCT(mosaic,.02);
    xExtra=1.2;  % extra abscissa in histogram
    cmExtra=.1;
    yExtra=1.1;  % extra headroom in histogram
    hMin=10;  % minimum log bin is 10 pixels
    ny=256;   % no. of y pixels in histogram
    minLx=.02;  % Min distance between limit lines
    
    fc=.1;
    colorMode=1;
    
    figure(1);
    set(gcf,'menubar','none');
    clf;
    sz=size(mosaic);
    nBin=2*round(max(sz)/2048);  % Binning factor, must be even
    nBin=max(nBin,1);
    nw=floor(sz/nBin);  % working image size
    smallMosaic=Crop(BinImage(mosaic,nBin),max(nw));  % make it square
    [nsx,nsy]=size(smallMosaic);
    
    % Create 256-entry color maps
    nj=size(jet,1);
    step=2*nj/256;
    cmShift=floor(cmExtra*2*nj);
    jets=circshift(jet,[-cmShift,0]);
    jets(nj-cmShift+1:nj,:)=repmat(jets(nj-cmShift,:),cmShift,1);  % saturate at top
    jeti=interp1(jets,step:step:nj,'linear',0);  % pad the lower end with zeros.
    jet0=[repmat(jeti(1,:),nj/2,1);jeti];
    
    gray0=repmat((0:255)'/256,1,3);
    % jet0=jet;
    % jet0(1,:)=0;  % lowest value is black
    % gray0=gray;
    
    cmSize=size(jet0,1);
    colormap(gray0);
    % colormap(gray);
    
    % Initialize the histogram
    mxLimit=Percentile(GaussFilt(smallMosaic,fc),.99999);  % Max intensity in image
    mxIntens=xExtra*mxLimit;  % Max intensity in histogram
    nBins=cmSize;   % no. of bins in intensity histogram = size of colormap
    sclb=nBins/mxIntens;  % scaling from image intensity to bin
    sMosaic=(GaussFilt(smallMosaic,fc)*sclb);  % image now has bin nos. as intensity.
    
    % Define two axes in the figure
    set(gcf,'name',fname);
    h2Height=.12;  % size of histogram port
    % Main image port
    h1=axes('position',[0.02 h2Height+.03 0.96 .96-h2Height]);
    imac(sMosaic);
    axis off
    % Histo port
    h2=axes('position',[0.02 0.02 0.96 h2Height]);
    
    bins=(0:nBins-1);  % Bin values are normalized to mxIntens
    h=hist(sMosaic(:),bins);
    hMax=max(h(nBins/2:end));  % max count we'll plot
    yExtent=[log10(hMin) yExtra*log10(hMax)];  % Histo vertical limits
    yh=ny*log10(max(hMin,h)/hMin)/diff(yExtent); % Histo plot values
    histImage=ColorHisto(yh,ny);
    xs=bins/nBins*xExtra;  % scale goes from 0 to just over 1 (xExtra)
    % imac(xs,yExtent,histImage);
    
    
    lx=[0.9 1.0];
    vx=0;
    vx0=0;
    b=1;
    cx=lx(1);
    cy=0;
    disp(' ');
    disp('press f to enter filter frequency (typ. .05 to .2)');
    disp(' o to open another file; q to quit.');
    
    while (lower(b)~='q') && lower(b)~='o' % q = quit; o=open
        %         disp([cx cy b]);
        switch b
            case 1
                if cx<2 && cy<10  % clck on histogram
                    dists=abs(lx-cx);
                    [d,j]=min(dists);
                    lx(j)=cx;
                    if lx(2)-lx(1)<minLx  % enforce a minimum distance
                        lx(3-j)=lx(j)+minLx*sign(1.5-j);
                    end;
                    %                     Update the colormap
                else  % clicked on image
                    ix=min(nsx,max(round(cx),1));
                    iy=min(nsy,max(round(cy),1));
                    vx=sMosaic(ix,iy)/cmSize*xExtra;
                    vx0=smallMosaic(ix,iy);
                    
                end;
                
                px=max(1,min(nBins,round(lx*nBins/xExtra)));  % convert x coord to bin number
                cMap=gray;
                cBins=px(1)+1:px(2);
                cMap(cBins,:)=0.8*jet0(cBins,:)+0.2*cMap(cBins,:);  % replace grayscale with spectrum here.
                colormap(cMap);
                %                 cmap(offs=cmSize*lx(1);
                %                 scl=1/diff(lx);
            case 'f'
                fc=MyInput('Filter freq ',fc);
                tMosaic=GaussFilt(smallMosaic,fc);
                mxLimit=Percentile(tMosaic,.99999);  % Max intensity in image
                mxIntens=xExtra*mxLimit;  % Max intensity in histogram
                sclb=nBins/mxIntens;  % scaling from image intensity to bin
                sMosaic=(tMosaic*sclb);  % image now has bin nos. as intensity.
                axes(h1);
                imac(sMosaic);
                h=hist(sMosaic(:),bins);
                hMax=max(h(nBins/2:end));  % max count we'll plot
                yExtent=[log10(hMin) yExtra*log10(hMax)];  % Histo vertical limits
                yh=ny*log10(max(hMin,h)/hMin)/diff(yExtent); % Histo plot values
                histImage=ColorHisto(yh,ny);
                
        end; % switch
        
        axes(h2);
        imac(xs,yExtent,histImage);
        hold on;
        for i=1:2
            plot([1 1]*lx(i),yExtent,'w-');
        end;
        plot([1 1]*vx,yExtent,'w--');
        hold off;
        if vx0>0
            text(double(vx),yExtent(2),['\color{white}' num2str(round(vx0))],...
                'VerticalAlignment','Top');
        end;
        %                 subplot(1,2,1);
        %                 imac((sMosaic-offs)*scl);
        
        [cx,cy,b] = Myginput(1,'square',1);
    end;
    ok=lower(b)=='o';
end;

