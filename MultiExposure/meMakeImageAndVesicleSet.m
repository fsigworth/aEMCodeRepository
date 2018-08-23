% meMakeImageAndVesicleSet.m
% phase-flip images and make a set of vesicle models

ds=2;
nZeros=inf;
shiftMax=10;  % max pixels shifted


% Ask for a set of original mi files.
[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
[rootPath infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;

cd(rootPath);
nfiles=numel(fname);
figure(1);
clf;
SetGrayscale;

for findex=1:nfiles
    % Get the mi structure
    load([infoPath fname{findex}]);
    
    mi0=mi;  % keep the old one
    %     update the weights
    mi.weights=weights;
    imgIndices=find(mi.weights>0);
    
    % Read the images
    nim=numel(mi.imageFilenames);
    fnames=cell(1,nim);
    disp('Processing images:');
    for j=1:nim
        fnames{j}=[mi.imagePath mi.imageFilenames{j}];
    end;
    
    m=meReadImages(fnames,cpe,0,mi.pixA,weights);
    m=mePreWhiten(m,modelSpectrum)/mi.doses(1);
    [nx ny nim]=size(m);
    n=[nx ny]/ds;

% Get individual images and CTFs, downsampled
    disp('Transforming images');

    [effctf, mc, mts, coeffz, dctfs]=meCombineImages(m,mi,ds,0,inf);
%     mts is the set of transformed images.).*
    %     Compute phase-flipped individual images.
    mSet=single(zeros([n nim]));  % Image set
    vSet=single(zeros([n nim]));  % Vesicle model set
    vRaw=single(zeros([n nim]));  % Vesicle model set
    
    shiftx=zeros(numel(mi.vesicle.x),numel(weights));
    shifty=zeros(numel(mi.vesicle.x),numel(weights));
    mis=mi;  % will be an array of structs.
%%    
    disp('Making phase-flipped images and vesicle models');
    for i=1:nim
        mis(i)=mi;
%         Zero out all the weights except for the i'th one.
        mis(i).weights=0*mi.weights;
        mis(i).weights(i)=1;
        indivCTF=dctfs(:,:,i);
%         Phase-flip the image
        mSet(:,:,i)=real(ifftn(ifftshift(-sign(indivCTF)).*fftn(mts(:,:,i))));
%         Shift the individual vesicle coordinates
        shiftx=mi.vesicle.shiftX(:,i);
        shifty=mi.vesicle.shiftY(:,i);
        truncShiftx=max(min(shiftx,shiftMax),-shiftMax);
        truncShifty=max(min(shifty,shiftMax),-shiftMax);

        mis(i).vesicle.x=mis(i).vesicle.x(:)+shiftx;
        mis(i).vesicle.y=mis(i).vesicle.y(:)+shifty;
        mis(i).vesicle.shiftx(:,i)=0;  % remove the shifts
        mis(i).vesicle.shifty(:,i)=0;
        
        v=1*meMakeModelVesicles(mis(i),n,0,0);  % no ctf
        vSet(:,:,i)=-real(ifftn(fftn(v).*(ifftshift(dctfs(:,:,i)))));
        vRaw(:,:,i)=v;
        
        subplot(221);
        imacs(v);
        title(fname{findex},'interpreter','none');
        subplot(222);
        imacs(vSet(:,:,i));
        title(i);
        subplot(223);
        imacs(BinImage(mSet(:,:,i),2));
        subplot(224);
        imacs(BinImage(mSet(:,:,i)-vSet(:,:,i),2));
        drawnow;
        
    end;
    freqs=ifftshift(RadiusNorm(n)/(mi.pixA*ds));
    [coeffs, mergedCTF, dctfs1]=meComputeMergeCoeffs( freqs, mi.ctf, mi.doses, 1, mi.weights);
    mMerged=-sum(real(ifft2(fft2(mts).*coeffs)),3);
    vMerged=sum(real(ifft2(fft2(vRaw).*dctfs1.*coeffs)),3);
    outName=[mi.procPath mi.baseFilename 'mvSet.mat'];
    save(outName,'mSet','vSet','mMerged','vMerged');
    disp(outName);
    jOutName=[mi.procPath 'jpeg/' mi.baseFilename];
    for i=1:nim
        WriteJpeg(mSet(:,:,i),[jOutName 'm' num2str(i) '.jpg']);
        WriteJpeg(vSet(:,:,i),[jOutName 'v' num2str(i) '.jpg']);
    end;
    WriteJpeg(mMerged,[jOutName 'mMerged.jpg']);
    WriteJpeg(vMerged,[jOutName 'vMerged.jpg']);
end;


