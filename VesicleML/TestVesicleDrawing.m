% TestVesicleDrawing.m

infoDir='Info_C24-4/';
tiffPath='Tiff_C24-4/';
miNames=f2FindInfoFiles(infoDir);
ds=4; % we know this is the image downsampling factor.

nmi=numel(miNames);
figure(1);
clf;

for i=1:1 % index of the micrograph
    mi=ReadMiFile(miNames{i});
    nv=numel(mi.vesicle.x); % Check if there are any vesicles.
    if nv<1
        continue;
    end;

    %     Create the affine matrix to handle scaling and image padding
    shft=(mi.imageSize-mi.padImageSize)/2;
    M=[ds 0 shft(1); 0 ds shft(2); 0 0 1];
    iM=inv(M);
    %     Get the local coordinates on our downsampled image according to
    %     local=[x;y;1]*iM;
    %     local is [ ix; iy; 1 ];
    %     but because Matlab uses 1-based arrays, we have to add 1 to ix and iy
    %     if we use them to index the image.

    imageName=[tiffPath mi.baseFilename 'ms.tif'];
    m=rot90(imread(imageName),-1); % rotation to correponds to EM convention
    % draw vesicle curves
    globalCoords=[mi.vesicle.x mi.vesicle.y ones(nv,1)];
    localCoords=iM*globalCoords'+1;
    lx=localCoords(1,:);
    ly=localCoords(2,:);
    lr=mi.vesicle.r/ds;

    subplot(2,1,1);
    imags(m);
    hold on;
    for j=1:nv
        [xs,ys]=CircleLineSegments((lr(j,:)));
        plot(xs+lx(j),ys+ly(j),'b-','linewidth',1);
    end;
    hold off;
    title([num2str(i) ':  ' imageName],'interpreter','none');
    drawnow;
%% we'll set single pixels in a blank image.
    m2=zeros(size(m),'single');
    [nx,ny]=size(m);
    subplot(2,1,2);
    for j=1:nv
%         Ask for a distance of 0.5 between points.
        [xs,ys]=CircleLineSegments(lr(j,:),0.5);
        ixs=round(xs+lx(j));
        iys=round(ys+ly(j));
        badLocs=(ixs<1)|(ixs>nx)|(iys<1)|(iys>ny);
        ixs(badLocs)=[];
        iys(badLocs)=[];
        for k=1:numel(ixs)
            m2(ixs(k),iys(k))=1;
        end;
    end;
    imags(m2);
    drawnow;

end;
