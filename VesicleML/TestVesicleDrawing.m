% TestVesicleDrawing.m
% Given Info/*mi.txt metadata files, and .tif images, create the
% corresponding vesicle models and center-of-membrane traces for the found
% vesicles.

numRots=-1; % default
writeSubtractions=1;
writeTracks=1;
fileSuffix='ms.tif';

% Set up directories for various datasets

% % The 20211122/C24-4 "low quality" dataset
% infoDir='/Users/fred/EMWork/Yangyu/20211122_farnam/C24-4_part/Info_C24-4/';
% tiffPath='/Users/fred/EMWork/Yangyu/20211122_farnam/C24-4_part/Tiff_C24-4/';
% fileSuffix='ms.tif';

% % The original Dec17 "high quality" dataset, ds=4
cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Data for Chris/181217ds4');
infoDir='Info/';
fileSuffix='ms.tif';
tiffPath='Micrographs/';
% 
% %  Vesicle-finder only fits of simulated vesicles
% cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Data for Chris/Chris sims VF only')
% infoDir='Info/';
% tiffPath='Micrographs/';
% numRots=0;

% % Refined fits of simulated vesicles
% cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Data for Chris/Chris sims/');
% infoDir='Info/';
% tiffPath='Micrographs/';
% numRots=0;

miNames=f2FindInfoFiles(infoDir);

nmi=numel(miNames);
figure(1);
clf;

for i=1:nmi % index of the micrograph
    disp([num2str(i) ':  ' miNames{i}]);
    mi=ReadMiFile(miNames{i});
    nv=numel(mi.vesicle.x); % Check if there are any vesicles.
    if nv<1
        continue;
    end;
    if ~isfield(mi,'padImageSize')
        mi.padImageSize=NextNiceNumber(mi.imageSize);
    end;

    imageName=[tiffPath mi.baseFilename fileSuffix];
    if ~exist(imageName,'file')
        disp(['File not found: ' imageName]);
        continue;
%         return
    end;
    m=rot90(imread(imageName),numRots); % rotation to correponds to EM convention
    ds=round(mi.padImageSize(1)/size(m,1));

    %     Create the affine matrix to handle scaling and image padding
    shft=(mi.imageSize-mi.padImageSize)/2;
    M=[ds 0 shft(1); 0 ds shft(2); 0 0 1];
    iM=inv(M);
    %     Get the local coordinates on our downsampled image according to
    %     local=[x;y;1]*iM;
    %     local is [ ix; iy; 1 ];
    %     but because Matlab uses 1-based arrays, we have to add 1 to ix and iy
    %     if we use them to index the image.

    % draw vesicle curves
    globalCoords=[mi.vesicle.x mi.vesicle.y ones(nv,1)];
    localCoords=iM*globalCoords'+1;
    lx=localCoords(1,:);
    ly=localCoords(2,:);
    lr=mi.vesicle.r/ds;

    mysubplot(2,1,1);
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
    mysubplot(2,1,2);
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
    if writeTracks
        trackName=[tiffPath mi.baseFilename '_track.tif'];
        imwrite(uint8(255*rot90(m2,0)),trackName,'TIFF');
    end;
    drawnow;

    if writeSubtractions
        disp('Computing vesicle models...')
        mq=meMakeModelVesicles(mi,[nx ny]);
        imwrite((uint8(imscale(mq,256))),[tiffPath mi.baseFilename 'memb.tif'],'TIFF');
        disp('done.')
    end;


pause(0.1)
end;
