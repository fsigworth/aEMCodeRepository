% meTrackVesicleShifts
% Given mi files, find vesicles in each of the underlying micrographs, and
% write this information into mi.vesicle.xshift and yshift

shiftMax=10;
cpe=16;
ds=2;  % Downsampling of merged image
iCamera=1;  % F20 ultrascan
weights=[1 1 1];

modelSpectrum=CCDModelSpectrum2D(iCamera);  % handle DE-12 or CCD

% Have the user select some mi files: boilerplate
if ~exist('fname','var') || ~exist('doBatchProcessing','var') || ~doBatchProcessing
    [fname, pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
    cd(rootPath);
end;

nfiles=numel(fname);

%%
for findex=1:nfiles
    figure(1);
clf;
SetGrayscale;

    % Get the mi structure
    load([infoPath fname{findex}]);
%%    
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
    
    %  Do the fitting here
    
    [effctf, mc, mts, coeffz, dctfs]=meCombineImages(m,mi,ds,0,inf);
    
    %     Compute phase-flipped individual images.
    mtsf=mts;  % To receive flipped images.
    mis=mi;
    shiftx=zeros(numel(mi.vesicle.x),numel(weights));
    shifty=shiftx;
    for i=imgIndices
        mis(i)=mi;
        mis(i).weights=0*mi.weights;
        mis(i).weights(i)=1;
        indivCTF=dctfs(:,:,i);
        mtsf(:,:,i)=real(ifftn(ifftshift(-sign(indivCTF)).*fftn(mts(:,:,i))));
        
        mis(i)=rsRefineVesicleFitsSub(mis(i),mtsf(:,:,i),2);

        shiftx(:,i)=(mis(i).vesicle.x(:)-mi.vesicle.x(:));
        shifty(:,i)=(mis(i).vesicle.y(:)-mi.vesicle.y(:));
    end;
%     truncShiftx=max(min(shiftx,shiftMax),-shiftMax);
%     truncShifty=max(min(shifty,shiftMax),-shiftMax);
    mi.vesicle.shiftX=shiftx;
    mi.vesicle.shiftY=shifty;
    mi.vesicle.shiftOk=abs(shiftx)<shiftMax & abs(shifty)<shiftMax;
    outName=[infoPath fname{findex}];
    save(outName,'mi');
    disp([outName ' saved.']);

    %%
    figure(1)
    subplot(1,1,1);
    imacs(GaussFilt(-mc,.1));
    drawnow;
    
    figure(2);
    mePlotVesicleShifts(mi);
    
end;
