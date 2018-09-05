function VesicleFinder(fileList,mpars);
% non-gui version of Vesicle_finding_GUI
% e.g. fileList={'ami.txt';'bmi.txt'};
% e.g. pars.sav.vesicleAmps=[6e-4 2e-3 0];

if nargin<2
    mpars=struct;
    mpars.sav=struct;
end;

disp('VesicleFinder startup.');



% Replace any given pars with the default amPars
% % pars=SetOptionValues(pars,mpars);

%% Context variables that are stored in a file
amPars.n0=256;  % Automask parameters
amPars.thresh=.1;
amPars.var=0;
amPars.width=0;
amPars.edge=.5;
amPars.denseFilt=.0015;  % A^-1
amPars.dense=.5;

sav.automaskPars=amPars;  % saved variables, stored in "VFContext.mat"
sav.basePath='';
sav.fullInfoPath='';
sav.baseName='';
sav.vesicleAmps=[2e-3 5e-3 0];
sav.vesicleRadii=[100 300 10];
sav.filterFreqs=[0 .1];   % HP and LP, in A^-1
sav.contrastPars=zeros(2,1);
sav.membranePars=[1.6 60 6];
sav.beamPars=[0 0 100];  % X, Y and R
sav.black=-2;
sav.white=.5;
sav.automaskFixed=0;
sav.initTheVesicles=1;
sav.eraseOldPicks=1;

% State variables
h.sav=sav;

h.batchMode=1;
h.fileList=fileList;
disp([num2str(numel(h.fileList)) ' files to process.']);

h.fileIndex=0;

h.imageLoaded=false;
h.readOnly=false;
h.pixA=1;
h.origImage=single(0);  % Original merged image
h.rawImage=single(0);   % merged image rescaled to display size
h.filtImage=single(0);  % rawImage after low and high pass filtering
% h.rawVesImage=single(0);
h.goodVesImage=single(0);
h.badVesImage=single(0);
h.ctf=single(0);
h.filtVesImage=single(0);
h.ccVals=0;     % values of cc maxima
h.ccValsScaled=0;
h.ccRadii=0;    % est. radius corresponding to cc maximum.
h.displayMode=0;  % show the raw data
h.maxDisplayMode=3;
h.oldFilterFreqs=[0 0];
h.maskIndex=3;       % default uses existing masks.
h.vesModels=single(0);
h.ds=1;
h.miChanged=false;
h.mi=struct;
h.miOriginal=struct;
h.miChanged=0;
h.automaskOn=false;
h.varianceMap=single(0);
h.axes1ButtonDown=false;
h.markedVesicleIndex=0;
h.findInMask=1;
% h.eraseOldPicks=0;
h.borderFraction=256;  % relative size of merged-image border is 1/256 ------------
h.theImage=[];
h.manualMaskDiameter=0;
h.manualMaskActive=false;
h.manualMaskCoords=zeros(0,'single');
h.oldVesicleModel=[];
h.e1CtrValue=0;
h.useFirstExposure=0;  % Flag to ignore first exposure in doing masking.
h.doTrackMembranes=1;
h.findingStarted=0;

% if ~isfield(h,'automaskBeamOn')
h.automaskBeamOn=0;  % don't change this if alreay defined.
% end;

h.displaySize=960;
f1=figure(1);
clf;
fpos=get(f1,'position');
set(f1,'position',[fpos(1:2) h.displaySize+[320+70 30] ]);

h.axes1=axes(f1, 'units','pixels','position',[10 10 960 960]);
% imags(randn(960));
h.axes3=axes(f1, 'units','pixels','position', [960+60 100 320 320]);
% plot(randn(10,1),randn(10,1),'o');

%%%%% Allow mpars to overwrite h values, except for h.sav
h=SetOptionValues(h,rmfield(mpars,'sav'));

%%%%% read a context file to get the parameters h.sav
% Try first the current directory, then the code directory.
pa=fileparts(which('Vesicle_finding_GUI'));
localName='VFContext.mat';
rootName=[AddSlash(pa) localName];
if exist('VFContext.mat','file');  % try to find the local one.
    disp(['Loading ' localName]);
    sav=load(localName);
elseif exist(rootName,'file')
    disp(['Loading ' rootName]);
    sav=load(rootName);
    disp(['Saving a local copy of ' localName]);
    save(localName,'sav');  % store a local copy.
else
    error([localName ' or ' rootName ': file not found']);
end;

%%%%%% Allow mpars.sav to overwrite the parameters from the file.
psav=SetOptionValues(sav.sav,mpars.sav);
h.sav=psav;

%     %%%%%%
%     h.sav.initTheVesicles=1;

h.imageFileTypes={'m.mrc' 'mf.mrc' 'm.jpg' 'mf.jpg'};  % allowed types of
% images to load based on mi file.

% startup
hObject=[];

h.jpegDir='Jpeg/';
CheckAndMakeDir(h.jpegDir,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robo-fit loop

for fIndex=1:numel(fileList)
    [h,imageOk]=LoadFile(h,fileList{fIndex});
    if imageOk
        disp('Automatic vesicle finding:');
        %         delete the old automask
        h=NewAutomask(h,false);
        disp('Automask off');
        h.doTrackMembranes=0;
        h=DoFind(h);
        h=DoFindMore(h);
        %%         Make a new automask
        disp('Automask on.');
        h=InitAutomask(h);
        h=NewAutomask(h,true);
        h=ShowImage(h);
        
        h.doTrackMembranes=1;
        h=DoFind(h);
        h=DoFindMore(h);
        
%         h=NewAutomask(h);
        h=CloseFile(h);
        jpegName=[h.jpegDir h.mi.baseFilename 'vf.jpg'];
        disp(['Saving ' jpegName]);
        print('-djpeg',jpegName);
    else
        disp(['No image found; skipping ' fileList{fIndex}]);
    end;
end;
disp('VesicleFinder finished.');
return  % ---------------end of Main----------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function h=CloseFile(h) % store the results of operations
        if ~h.imageLoaded || ~h.miChanged
            disp('No changes to mi file.  Not written.');
            return
        end;
        
        if h.readOnly
            disp('Read-only mode.  No file written.');
            return
        end;
        
        % discard unused masks
        if isfield(h.mi,'mask')&&(numel(h.mi.mask)>h.maskIndex)
            h.mi.mask=h.mi.mask(1:h.maskIndex);
        end;
        
        if h.sav.eraseOldPicks
            disp('Erasing old particles.');
            h.mi.particle.picks=[];
            mi=h.mi;
        else
            disp('Merging the vesicle list.');
            mi=rsMergeVesicleList2(h.mi,h.miOriginal);
        end;
        mi=RemoveOverlappingVesicles(mi);
        if numel(h.oldVesicleModel)>0  % Restore the old vesicle model
            mi.vesicleModel=h.oldVesicleModel;
        end;
        mi.log{end+1,1}=['VesicleFinder ' TimeStamp];
        outName=[mi.infoPath mi.baseFilename 'mi.txt'];
        nameWritten=WriteMiFile(mi,outName);
        disp(['Saved: ' nameWritten]);
    end


    function [h,imageOk]=LoadFile(h,fileName)
        
        h=CloseFile(h);  % Save the previous file.
        h.miChanged=0;
        
        disp(' ');
        disp(['Reading ' fileName]);
        mi=ReadMiFile(fileName);  % loads mi
        nim=min(min(numel(mi.ctf),numel(mi.doses)),size(mi.frameSets,1));
        for i=1:nim
            disp([num2str(mi.ctf(i).defocus,3) 'um.  frames: ' num2str(mi.frameSets(i,:)) '  dose: ' num2str(mi.doses(i),3)]);
        end;
        h.sav.baseName=mi.baseFilename;
        
%         if h.sav.initTheVesicles
%             mi=ZeroOutVesicles(mi,h);
%         end;
%         
        [h, imageOk]=GetImageFile(mi,h);
        h.imageLoaded=imageOk;
        if ~imageOk
            return
        end;
        h=InitDisplay(h);
        drawnow;
    end


    function mis=ZeroOutVesicles(mis,h)
        % Initialize all the vesicle fields in the mi structure
        % Also initializes the membrane model
        vLipid=h.sav.membranePars(1);
        thk=h.sav.membranePars(2);
        rise=h.sav.membranePars(3);
        pixA=mis.pixA;
        % Create the model, which is sampled in units of the original pixel size.
        nm0=ceil(thk/(2*pixA))*2+3;  % array for vesicle model; 60A nominal
        mis.vesicleModel=fuzzymask(nm0,1,thk/pixA/2,rise/pixA)...
            *vLipid;  % units of V
        % Zero out the previous picks
        mis.vesicle.x=[];
        mis.vesicle.y=[];
        mis.vesicle.r=[];
        mis.vesicle.s=[];
        mis.vesicle.ok=[];
        mis.vesicle.shiftX=[];
        mis.vesicle.shiftY=[];
        mis.vesicle.shiftOk=[];
        mis.vesicle.af=[];
        mis.vesicle.refined=0;
        mis.vesicle.extraPeaks=[];
        mis.vesicle.extraSD=0;
        mis.vesicle.extraS=[];
    end


    function [h,ok]=GetImageFile(mi,h)
        % Having loaded an mi file, load the image and initialize variables.
        
        % First, update the vesicle info to the canonical nv x 4 logical matrix
        nv=numel(mi.vesicle.x);
        if ~(isfield(mi.vesicle,'ok')...
                && size(mi.vesicle.ok,1)==nv)
            mi.vesicle.ok=true(nv,1);
        end;
        if size(mi.vesicle.ok,2)<4 && nv>0
            mi.vesicle.ok(nv,4)=false;  % extend the array
        end;
        
        % Load the merged image
        imageBasename=[mi.procPath mi.baseFilename 'm.mrc'];
        % ok=0;
        %
        sufExts={'s.mrc' 'z.tif' '.mrc'};
        [fullImageName,ok]=CheckForAltImage(imageBasename,sufExts);  % valid filename?  Load it
        h.oldFilterFreqs=[0 0];
        h.baseName=imageBasename;
        
        if ok
            %     h.origImage=ReadEMFile(fullImageName)/mi.doses(1);  % normalize by dose.
            if strcmp(fullImageName(end-4:end),'z.tif')
                disp('Reading the compressed merged image');
            end;
            h.origImage=ReadEMFile(fullImageName);
            
            %     Pick up the first exposure for use in masking
            exp1name=mi.imageFilenames{1};  % get the low-defocus image
            [exp1name2,ok1]=CheckForImageOrZTiff([mi.imagePath exp1name]);
            if ok1 && h.useFirstExposure  % pick up the first exposure and make it the same size.
                if ~strcmp(exp1name2(end-4:end),'z.tif')
                    disp('Reading the first exposure');
                else
                    disp('Reading the compressed first exposure');
                end;
                e1Image=ReadEMFile(exp1name2);
            else
                if h.useFirstExposure
                    disp(['First exposure not found: ' exp1name]);
                    disp(' ...using merged image for an inferior global mask.');
                end;
                e1Image=h.origImage+5;  % offset it from zero.
            end;
            h.exp1Image=DownsampleGeneral(e1Image,h.displaySize/2);
            
            h.imageLoaded=true;
            h.findingStarted=false;
            h.rawImage=single(0);
            h.automaskOn=false;
            %     Clear all the processed images.
            h.filtImage=single(0);
            h.rawVesImage=single(0);
            h.filtVesImage=single(0);
            
            % initialize the mask
            t=mi.mergeMatrix;
            msk=meMakeMergedImageMask(mi.imageSize/4,t,mi.imageSize/(4*h.borderFraction));
            mi=meInsertMask(msk,mi,1);
            %     if h.automaskBeamOn  % we've left this on from last time; update values
            %         mi.mask(2).merge='AND';
            %         mi.mask(2).encoding='RIM';
            %         mi.mask(2).data=h.sav.beamPars/100;
            %     end;
            %     Point to the end of the stack
            h.maskIndex=numel(mi.mask);
            
            %   Initialize variables
            h.oldVesicleModel=mi.vesicleModel;
            h.miOriginal=mi;
            h.mi=mi;  % -------copy the mi structure here------
            h.ccVals=0;
            h.ccValsScaled=0;
            h.ccRadii=0;
            h.markedVesicleIndex=0;
            h.ifImageComp=zeros(256,'single');
            %     h.ifImageFlat=zeros(256,'single');
            %     if h.maskIndex<3 % no audomasking has been done
            %         set(h.togglebutton_Automask,'value',true);
            %     end;
            % else
            %     msgbox(['Can''t find the image ' imageBasename],'ok');
        end;
    end


% ----InitDisplay----
    function h=InitDisplay(h)
        if ~h.imageLoaded
            return
        end;
        % set(h.figure1,'pointer','watch');
        osize=size(h.origImage);
        %         Get the downsampling factor from the loaded image, and from the
        %         original micrograph.
        h.dsImage=(max(osize./h.displaySize)); % Force image to fit
        h.ds0=h.dsImage*h.mi.imageSize(1)/osize(1);
        h.pixA=h.ds0*h.mi.pixA; % scale of displayed image.
        
        %         Create the displayed image
        h.rawImage=DownsampleGeneral(h.origImage,osize(1)/h.dsImage);
        
        h=CreateE1Map(h);
        
        disp('Making model vesicles')
        if ~isfield(h.mi.vesicle,'ok')
            h.mi.vesicle.ok=false;  % make one entry.
        end;
        [nv, ne]=size(h.mi.vesicle.ok);
        if ne<4
            h.mi.vesicle.ok(1,4)=false;  % extend it.
        end;
        
        goodVes=all(h.mi.vesicle.ok(:,1:2),2); % vesicles in range
        badVes=(h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2)); % found, but not in range
        
        h.goodVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(goodVes));
        h.badVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(badVes));
        disp('...done');
        h.ctf=meGetEffectiveCTF(h.mi,size(h.rawImage),h.ds0);
        h=UpdateDisplayFiltering(h);
        % h=UpdateAutomaskBeam(h);  % Also calls ShowImage.
    end


% UpdateDisplayFiltering
    function h=UpdateDisplayFiltering(h,doRawImage)
        if nargin<2
            doRawImage=1;
        end;
        h.filtImage=h.rawImage;
        h.filtVesImage=h.goodVesImage+h.badVesImage;
        
        fcs=h.sav.filterFreqs*h.pixA;
        if fcs(1)  % highpass freq not zero
            if doRawImage
                h.filtImage=GaussHP(h.filtImage,fcs(1));
            end;
            h.filtVesImage=GaussHP(h.filtVesImage,fcs(1));
        end;
        if fcs(2)  % lowpass freq not zero
            if doRawImage
                h.filtImage=GaussFilt(h.filtImage,fcs(2));
            end;
            h.filtVesImage=GaussFilt(h.filtVesImage,fcs(2));
        end;
        h.oldFilterFreqs=h.sav.filterFreqs;
    end



% ----ShowImage----
    function h=ShowImage(h)
        if h.imageLoaded
            
            % Draw the scatterplot of existing vesicles
            if isfield(h.mi.vesicle,'ok')
                [nv, ne]=size(h.mi.vesicle.ok);
                if nv>0 && ne>3 && numel(h.mi.vesicle.r)>0
                    goodVes=all(h.mi.vesicle.ok(:,1:2),2); % vesicles in range
                    badVes=(h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2)); % found, but not in range
                    plot(h.axes3,...
                        h.mi.vesicle.r(goodVes,1)*h.mi.pixA,h.mi.vesicle.s(goodVes,1),'b.',...
                        h.mi.vesicle.r(badVes,1)*h.mi.pixA,h.mi.vesicle.s(badVes,1),'m.',...
                        'markersize',15);
                    if h.markedVesicleIndex>0
                        set(h.axes3,'nextplot','add');
                        plot(h.axes3,h.mi.vesicle.r(h.markedVesicleIndex,1)*h.mi.pixA,...
                            h.mi.vesicle.s(h.markedVesicleIndex,1),'ks','markersize',12);
                        set(h.axes3,'nextplot','replace');
                    end;
                    title(h.axes3,[num2str(sum(goodVes)) ' vesicles in range']);
                else
                    title(h.axes3,' '); % clear the title
                    goodVes=false;
                    badVes=false;
                end;
            end;
            nv=numel(h.mi.vesicle.x);
            
            % if numel(h.filtImage)>0  % an image has been loaded
            msk=rot90(meGetMask(h.mi,size(h.filtImage),1:h.maskIndex));
            imData=h.filtImage;
            n=size(h.filtImage);
            showMask=1;
            showAmps=0;
            showCircles=0;
            showGhosts=0;
            switch h.displayMode
                case 0
                    imData=h.filtImage;
                    showCircles=1;
                    showAmps=1;
                case 1
                    imData=h.filtImage-h.filtVesImage;
                case 2
                    imData=h.filtImage-h.filtVesImage;
                    showGhosts=1;
                case 3
                    imData=h.filtImage-h.filtVesImage;
                    showGhosts=1;
                    showAmps=1;
                case 4
                    if numel(h.ccValsScaled)>1
                        imData=h.ccValsScaled(1:n(1),1:n(2));
                        showGhosts=1;
                    else
                        return
                    end;
                case 5
                    imData=Downsample(h.ifImageComp,size(h.filtImage));
                    showGhosts=1;
            end;
            %     theImage =  repmat(rot90(imscale(imData,256,1e-3)),[1 1 3]);
            midValue=h.e1CtrValue/(h.mi.doses(1)*h.mi.cpe)-1;
            %     theImage =  repmat(rot90(256*(imData-midValue-h.sav.black)/(h.sav.white-h.sav.black)),[1 1 3]);
            % We just autoscale the image.
            theImage =  repmat(rot90(imscale(imData,256,.0001)),[1 1 3]);
            
            nx=size(h.rawImage,1);
            ny=size(h.rawImage,2);
            
            if showGhosts
                ghostColor=[.7 .7 1];
                ghostColorBad=[1 .5 .35];
                dotW=2;
                
                xCtrs=max(dotW+1,min(nx-dotW,round(h.mi.vesicle.x/h.ds0+1)));
                yCtrs=max(dotW+1,min(ny-dotW,round(h.mi.vesicle.y/h.ds0+1)));
                if any(goodVes)
                    gves=imscale(max(-h.goodVesImage,-1e-9),1);
                    for i=find(goodVes)'
                        gves(xCtrs(i)-dotW:xCtrs(i)+dotW,yCtrs(i)-dotW:yCtrs(i)+dotW)=1.5;
                    end;
                else
                    gves=0;
                end;
                if any(badVes)
                    bves=imscale(max(-h.badVesImage,-1e-9),1);
                    for i=find(badVes)'
                        bves(xCtrs(i)-dotW:xCtrs(i)+dotW,yCtrs(i)-dotW:yCtrs(i)+dotW)=1.5;
                    end;
                    
                    %             hp1=plot(h.mi.vesicle.x(goodVes)/h.ds0+1,...
                    %             ny-(h.mi.vesicle.y(goodVes)/h.ds0),'b.','markersize',10);
                    
                else
                    bves=0;
                end;
                
                color=(1-ghostColor); % color to subtract for membrane
                for i=1:3
                    theImage(:,:,i)=theImage(:,:,i).*(1-rot90(gves)*color(i));
                end;
                color=(1-ghostColorBad); % color for bad vesicle
                for i=1:3
                    theImage(:,:,i)=theImage(:,:,i).*(1-rot90(bves)*color(i));
                end;
            end;
            if showMask
                maskColor=[1 .8 .85];
                
                color=(1-maskColor);
                for i=1:3
                    theImage(:,:,i)=theImage(:,:,i).*(1-(1-msk)*color(i));
                end;
            end;
            %     --------draw the image------
            axes(h.axes1);
            cla reset
            ih = imshow(uint8(theImage),'InitialMagnification',100,'parent',h.axes1);
            set(h.axes1,'xLimMode','manual');
            set(h.axes1,'yLimMode','manual');
            title(h.mi.baseFilename,'interpreter','none');
            %     if showGhosts % draw the center points of the vesicles
            %         %         goodCenterColor=[.5 .5 1];
            %         %         badCenterColor=[1 .5 1];
            %         if isfield(h.mi.vesicle,'ok') && size(h.mi.vesicle.ok,2)>3
            %             badVes=h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2);
            %             goodVes=all(h.mi.vesicle.ok(:,1:2),2);
            %         else
            %             badVes=false(numel(h.mi.vesicle.x),1);
            %             goodVes=true(numel(h.mi.vesicle.x),1);
            %         end;
            %         hold on;
            %         ny=size(h.rawImage,2);
            %         hp1=plot(h.mi.vesicle.x(goodVes)/h.ds0+1,...
            %             ny-(h.mi.vesicle.y(goodVes)/h.ds0),'b.','markersize',10);
            %         hp2=plot(h.mi.vesicle.x(badVes)/h.ds0+1,...
            %             ny-(h.mi.vesicle.y(badVes)/h.ds0),'r.','markersize',10);
            %         hold off;
            %         set(hp1,'HitTest','off');
            %         set(hp2,'HitTest','off');
            %     end;
            if showAmps
                for i=1:nv
                    x=double(h.mi.vesicle.x(i)/h.ds0+1);
                    y=double(ny-h.mi.vesicle.y(i)/h.ds0+1);
                    amp=h.mi.vesicle.s(i,1)*1000;
                    text(x,y,num2str(amp,2),'color',[.8 .8 0],'fontsize',12,'PickableParts','none');
                end;
            end;
            
            if showCircles
                hold on
                for i=1:nv
                    r1=h.mi.vesicle.r(i,:)/h.ds0;
                    [x,y]=CircleLineSegments(h.mi.vesicle.r(i,:)/h.ds0,min(10,100/r1(1)));
                    x=double(x+h.mi.vesicle.x(i)/h.ds0+1);
                    y=double(y+ny-h.mi.vesicle.y(i)/h.ds0+1);
                    if goodVes(i)
                        plot(x,y,'b-','HitTest','off');
                    elseif badVes(i)
                        plot(x,y,'r-','HitTest','off');
                    end;
                end;
                hold off
            end;  % if show
            drawnow;
        end % if h.imageLoaded;
    end


%     function [h, ok]=LoadAnotherMiFile(h,offset)
%         ok=0;
%         if ~h.imageLoaded
%             return
%         end;
%         CloseFile(h);  % Batch mode: save the previous one.
%         h.miChanged=0;
%         h.fileIndex=h.fileIndex+offset;
%         infoPath='';
%         if h.fileIndex> numel(h.fileList)
%             return
%         end;
%         fileName=h.fileList{h.fileIndex};
%         disp(['Reading ' fileName]);
%         [mi,nameRead]=ReadMiFile([infoPath fileName]);
%         
%         nim=min(min(numel(mi.ctf),numel(mi.doses)),size(mi.frameSets,1));
%         for i=1:nim
%             disp([num2str(mi.ctf(i).defocus,3) 'um.  frames: ' num2str(mi.frameSets(i,:)) '  dose: ' num2str(mi.doses(i),3)]);
%         end;
%         
%         [pa,nm,ex]=fileparts(nameRead);
%         fileName=[nm ex];
%         %     set(h.FileNameText,'String',fileName);
%         
%         h.sav.baseName=nm(1:numel(nm)-2);  % delete the 'mi'
%         mi.basePath=h.sav.basePath;
%         
%         [h ok]=GetImageFile(mi,h); % copies mi into h.mi
%         if ~ok
%             return
%         end;
%         
%         % Update the mask
%         
%         % if ~isfield(mi,'mask') || ~isfield(mi.mask,'merge') || numel(mi.mask.merge)<1
%         %    Old merged data: compute the base mask and insert it.
%         % disp('Computing merge mask...');
%         t=h.mi.mergeMatrix;
%         msk=meMakeMergedImageMask(h.mi.imageSize/4,t,h.mi.imageSize/(4*h.borderFraction));
%         h.mi=meInsertMask(msk,h.mi,1);
%         % end;
%         %     Point to the end of the stack
%         h.maskIndex=numel(h.mi.mask);
%         h=InitDisplay(h);
%         
%     end
% 
% 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------Automask------------

% --- Executes on button press in togglebutton_Automask.
    function togglebutton_Automask_Callback(hObject, eventdata, h)
        h.automaskOn=get(h.togglebutton_Automask,'value') && h.imageLoaded;
        if h.automaskOn  % First, compute the inverse-filtered images
            h=InitAutomask(h);
            h=NewAutomask(h);
        else
            h=NewAutomask(h,false);
        end;
        guidata(hObject, h);
    end

    function togglebutton_Automask_ButtonDownFcn(hObject,eventdata, h)
        h.automaskOn=get(h.togglebutton_Automask,'value') && h.imageLoaded;
        if h.automaskOn  % First, compute the inverse-filtered images
            h=InitAutomask(h);
            h=NewAutomask(h);
        else
            h=NewAutomask(h,false);
        end;
        guidata(hObject, h);
    end

% --- Executes on button press in fixedCheckbox.
    function fixedCheckbox_Callback(hObject, eventdata, h)
        h.sav.automaskFixed=get(hObject,'value');
        h=CreateE1Map(h);
        if h.automaskOn
            h=NewAutomask(h);
        end;
        guidata(hObject,h);
    end

% hObject    handle to fixedCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixedCheckbox

    function h=CreateE1Map(h)
        %   Create the global "dense" mask
        if numel(h.exp1Image)>1 && ~any(isnan(h.exp1Image(:))) % Something there
            %     filter it, make a histogram, set the threshold
            exp1f=GaussFiltDCT(h.exp1Image,h.sav.automaskPars.denseFilt*h.pixA*2);
            med=abs(median(exp1f));  % make sure it's positive, should be.
            bins=.75*med:med/500:1.25*med;
            e1Hist=hist(exp1f(:),bins);
            [e1Val,e1Mode]=max(e1Hist);
            e1Upper=find(e1Hist>e1Val/2,1,'last');
            %     sigma5=5*(e1Upper-e1Mode+.5);
            e1Ctr=bins(e1Mode);
            e1Ctr=bins(e1Upper); %%%%%%%%
            if h.sav.automaskFixed
                e1Ctr=h.e1CtrValue;
            else
                h.e1CtrValue=e1Ctr;
            end;
            
            %     sigma1=std(exp1f(abs(exp1f-e1Ctr)<sigma5));
            %     e1mx=e1Ctr+sigma1;
            %     e1mn=e1Ctr-10*sigma1;
            %     h.e1Map=(exp1f-e1mn)/(e1mx-e1mn);
            h.e1Map=(10*(exp1f/e1Ctr-.9)); % expand the top 10%
        else
            h.e1Map=ones(size(h.exp1Image),'single');  % all ones
        end;
        % Also create the varMap.
        f1=(h.pixA*2)/20;  % lowpass cutoff of bp filter for variance
        f2=(h.pixA*2)/60;  % highpass cutoff
        varFilt=(h.pixA*2)/400; % filtering of the variance
        workImage=DownsampleGeneral(h.rawImage,h.displaySize/2);
        varMap=GaussFiltDCT((SharpHP(SharpFilt(workImage,f1),f2)).^2,varFilt);
        
        mxVal=max(varMap(:));  % make sure it's positive, should be.
        bins=0:mxVal/100:mxVal;
        varHist=hist(varMap(:),bins);
        [~,modeIndex]=max(varHist);
        varMode=bins(modeIndex);
        h.varMap=varMap/(10*varMode);
    end

    function h=InitAutomask(h)
        h.maskIndex=max(h.maskIndex,3);  % Show this mask
        m=h.rawImage-h.goodVesImage-h.badVesImage;
        if all(size(m)>1) % an actual image
            m=Downsample(m,size(m)/2);
            
            h.ifImage=meCTFInverseFilter(m,h.mi,1,0,0);  % totally inverse filtered
            % h.ifImageFlat=GaussFilt(meCTFInverseFilter(m,h.mi,1,0,.0005),.05);
            h.ifImageComp=GaussFilt(meCTFInverseFilter(m,h.mi,1,.0005,.0005),.1);
        else
            disp(['*** size of m is ' num2str(size(m))]);
        end;
    end

    function h=NewAutomask(h,active)
        if nargin<2
            active=true;
        end;
        if ~active
            h.mi=meInsertMask(true(round(size(h.rawImage)/2)),h.mi,3,'AND');
            h.miChanged=1;
            h=ShowImage(h);
            return
        end;
        % if h.automaskOn
        %     thr1=h.sav.automaskPars.thresh*(h.automaskLims(2)-h.automaskLims(1))...
        %         +h.automaskLims(1);
        thr2=2*h.sav.automaskPars.width-2;
        thr3=h.sav.automaskPars.dense;
        mskLocal=h.ifImageComp>thr2;
        %     mskGlobal= h.e1Map>thr3;
        mskGlobal= (h.e1Map>thr3) & (h.varMap<(1-h.sav.automaskPars.var));
        fc1=exp(-h.sav.automaskPars.edge*3);
        fc2=fc1/5;
        if fc1<.9
            mskLocal=GaussFilt((GaussFilt(mskLocal,fc1)>.1),fc2)>.9;
        end;
        fc1=exp(-h.sav.automaskPars.thresh*3)*.1;
        mskGlobal=(GaussFiltDCT(mskGlobal,fc1)>.99);
        
        % The automask always goes into position 3.
        h.mi=meInsertMask(mskLocal & mskGlobal,h.mi,3,'AND');
        h.miChanged=1;
        h=ShowImage(h);
        % end;
    end


% --------------------------finding---------------------------

    function h=DoFind(h)
        if h.sav.initTheVesicles && ~h.findingStarted
            h.mi=ZeroOutVesicles(h.mi,h);
            disp('Erasing old vesicles.');
        end;
        h.findingStarted=true;
        
        if ~h.imageLoaded
            disp('No image.');
            return
        end;
        
        % membrane model
        vLipid=h.sav.membranePars(1);
        thk=h.sav.membranePars(2);
        rise=h.sav.membranePars(3);
        pixA=h.mi.pixA;
        
        % Create the model, which is sampled in units of the original pixel size.
        nm0=ceil(thk/(2*pixA))*2+3;  % array for vesicle model; 60A nominal
        if ~isfield(h.mi,'vesicleModel') || numel(h.mi.vesicleModel)<3
            h.mi.vesicleModel=fuzzymask(nm0,1,thk/pixA/2,rise/pixA)...
                *vLipid;  % units of V
        end;
        rPars=h.sav.vesicleRadii;
        % Zero out the previous picks
        mi1=h.mi;
        mi1.vesicle.x=[];
        mi1.vesicle.y=[];
        mi1.vesicle.r=[];
        mi1.vesicle.s=[];
        mi1.vesicle.ok=[];
        if isfield(mi1,'mask')
            mi1.mask=mi1.mask(1:h.maskIndex); % collapse the mask array
        end;
        % Initialize the vesicle finder
        mi1=rsFindVesicles3(h.rawImage, mi1, rPars, h.findInMask);
        %%
        minAmp=h.sav.vesicleAmps(1);
        maxAmp=h.sav.vesicleAmps(2);
        
        mins=1;
        nVesOld=0;
        % Loop through finding groups of 50 vesicles.
        while mins>minAmp
            [mi1, t]=rsFindVesicles3('next',50,h.sav.vesicleAmps);
            axes(h.axes1);
            ih = imshow(rot90(uint8(imscale(t.ms-t.umodel,256,1e-3))),...
                'InitialMagnification',100*size(t.ms,1)/size(h.rawImage,1) );
            %     drawnow;
            mins=t.globalmax;
            nves=size(mi1.vesicle.s,1);
            %     set(hObject,'string',num2str(nves));
            %     Make the scatterplot of vesicle radii and amplitudes
            xs1=rPars(1);
            xs2=rPars(2);
            xs=[xs1    xs1    xs2    xs2    xs1];
            ys=[minAmp maxAmp maxAmp minAmp minAmp];
            if numel(mi1.vesicle.ok)<4  % not even one row present
                mi1.vesicle.ok=false(1,4); % make sure there are enough elements
            end;
            goodVes=all(mi1.vesicle.ok(1:nves,1:2),2);
            badVes=mi1.vesicle.ok(1:nves,1) & ~mi1.vesicle.ok(1:nves,2);
            %     plot(h.axes3,mi1.vesicle.r(1:nves)*mi1.pixA,mi1.vesicle.s(1:nves),...
            %         'b.', xs,ys,'r-','markersize',10);
            if numel(mi1.vesicle.r)>0  % something to draw
                plot(h.axes3,mi1.vesicle.r(goodVes,1)*mi1.pixA,mi1.vesicle.s(goodVes,1),'b.',...
                    mi1.vesicle.r(badVes,1)*mi1.pixA,mi1.vesicle.s(badVes,1),'m.',...
                    xs,ys,'k-','markersize',10);
            else
                cla(h.axes3);
            end;
            title(nves);
            %     Exit the loop when we can't find any more vesicles
            if nves<=nVesOld
                break
            end;
            nVesOld=nves;
        end;
        h.ccVals=t.ccmx;
        h.ccValsScaled=t.ccmxScaled;
        h.ccRadii=(t.fitmin+(t.ccmi-1)*t.rstep)*t.ds;  % radius in orig pixels
        
        mi1.vesicle.shiftX=[];
        mi1.vesicle.shiftY=[];
        mi1.vesicle.shiftOk=[];
        mi1.vesicle.af=[];
        mi1.vesicle.refined=0;
        mi1.vesicle.extraPeaks=[];
        mi1.vesicle.extraSD=0;
        mi1.vesicle.extraS=[];
        h.mi=mi1;
        h.miChanged=1;
        h.goodVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(goodVes));
        h.badVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(badVes));
        h.rawVesImage=h.goodVesImage+h.badVesImage;
        % h.rawVesImage=meMakeModelVesicles(h.mi,size(h.rawImage));
        h=UpdateDisplayFiltering(h);
        % h.displayMode=0;  % mark the vesicles
        h.markedVesicleIndex=0;
        ShowImage(h);
        drawnow;
        
        if h.doTrackMembranes
            
            % Tune up the vesicle fits
            disp('Tracking vesicle membranes');
            
            h.mi=TrackVesicleMembrane(h);
            
            % display everything again
            disp('Computing vesicle models');
            goodVes=all(h.mi.vesicle.ok(:,1:2),2); % vesicles in range
            badVes=(h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2)); % found, but not in range
            h.goodVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(goodVes));
            h.badVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(badVes));
            h.rawVesImage=h.goodVesImage+h.badVesImage;
            % h.rawVesImage=meMakeModelVesicles(h.mi,size(h.rawImage));
            h=UpdateDisplayFiltering(h);
            ShowImage(h);
            drawnow;
            
        end;
        disp('Done.');
    end


% --- Executes on button press in pushbutton_FindMore.-----------------
% Searches the subtracted micrograph for more vesicles
    function h=DoFindMore(h)
        if ~h.imageLoaded
            return
        end;
        
        % membrane model
        vLipid=h.sav.membranePars(1);
        thk=h.sav.membranePars(2);
        rise=h.sav.membranePars(3);
        pixA=h.mi.pixA;
        % Create the model, which is sampled in units of the original pixel size.
        nm0=ceil(thk/(2*pixA))*2+3;  % array for vesicle model; 60A nominal
        if ~isfield(h.mi,'vesicleModel') || numel(h.mi.vesicleModel)<3
            h.mi.vesicleModel=fuzzymask(nm0,1,thk/pixA/2,rise/pixA)...
                *vLipid;  % units of V
        end;
        rPars=h.sav.vesicleRadii;
        % in this case we use the raw vesicle image for finding.
        mi1=rsFindVesicles3(h.rawImage-h.rawVesImage, h.mi, rPars, h.findInMask);
        prevNFound=numel(mi1.vesicle.x);
        %%
        minAmp=h.sav.vesicleAmps(1);
        maxAmp=h.sav.vesicleAmps(2);
        mins=1;
        nVesOld=0;
        % Loop through finding groups of 50 vesicles.
        while mins>minAmp
            [mi1, t]=rsFindVesicles3('next',50,h.sav.vesicleAmps);
            axes(h.axes1);
            ih = imshow(rot90(uint8(imscale(t.ms-t.umodel,256,1e-3))),...
                'InitialMagnification',100*size(t.ms,1)/size(h.rawImage,1) );
            %     drawnow;
            mins=t.globalmax;
            nves=size(mi1.vesicle.s,1);
            %     set(hObject,'string',num2str(nves));
            %     Make the scatterplot of vesicle radii and amplitudes
            xs1=rPars(1);
            xs2=rPars(2);
            xs=[xs1    xs1    xs2    xs2    xs1];
            ys=[minAmp maxAmp maxAmp minAmp minAmp];
            if numel(mi1.vesicle.ok)<4
                mi1.vesicle.ok=false(1,4); % make sure there are enough elements
            end;
            goodVes=all(mi1.vesicle.ok(1:nves,1:2),2);
            badVes=mi1.vesicle.ok(1:nves,1) & ~mi1.vesicle.ok(1:nves,2);
            %     plot(h.axes3,mi1.vesicle.r(1:nves)*mi1.pixA,mi1.vesicle.s(1:nves),...
            %         'b.', xs,ys,'r-','markersize',10);
            if numel(mi1.vesicle.r)>0  % something to draw
                plot(h.axes3,mi1.vesicle.r(goodVes,1)*mi1.pixA,mi1.vesicle.s(goodVes,1),'b.',...
                    mi1.vesicle.r(badVes,1)*mi1.pixA,mi1.vesicle.s(badVes,1),'m.',...
                    xs,ys,'k-','markersize',15);
            else
                cla(h.axes3);
            end;
            title(nves);
            %     Exit the loop when we can't find any more vesicles
            if nves<=nVesOld
                break
            end;
            nVesOld=nves;
        end;
        h.ccVals=t.ccmx+h.ccVals;
        h.ccValsScaled=t.ccmxScaled+h.ccValsScaled;
        h.ccRadii=(t.fitmin+(t.ccmi-1)*t.rstep)*t.ds;  % radius in orig pixels
        
        mi1.vesicle.shiftX=[];
        mi1.vesicle.shiftY=[];
        mi1.vesicle.shiftOk=[];
        mi1.vesicle.af=[];
        mi1.vesicle.refined=0;
        mi1.vesicle.extraPeaks=[];
        mi1.vesicle.extraSD=0;
        mi1.vesicle.extraS=[];
        
        % Additional vesicles are marked 'false' in the 3rd-4th column.
        totalNFound=numel(mi1.vesicle.x);
        if totalNFound>prevNFound
            disp([num2str(totalNFound-prevNFound) ' additional vesicles found...']);
            mi1.vesicle.ok(prevNFound+1:totalNFound,3:4)=false;
        end;
        
        h.mi=mi1;
        h.miChanged=1;
        h.goodVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(goodVes));
        h.badVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(badVes));
        h.rawVesImage=h.goodVesImage+h.badVesImage;
        % h.rawVesImage=meMakeModelVesicles(h.mi,size(h.rawImage));
        h=UpdateDisplayFiltering(h);
        % h.displayMode=0;  % subtract and mark the vesicles
        h.markedVesicleIndex=0;
        ShowImage(h);
        disp('Done.');
    end


%
end
