function BatchRSPicker(miNames)
% Batch version of SimpleRSPicker
% F. Sigworth, Jan '23
%  See rspLoadPicksFromMi.m to see assignment of ptrs
%  See rspNewBox to see assignment of flag values
si=struct;

version=102;
disOk=0;


% Retrieve parameters from a file in the program directory
pa=fileparts(which('SimpleRSPicker'));
datName=[AddSlash(pa) 'SimpleRSPickerDat.mat'];
if exist(datName,'file')>0
    dis=load(datName);
    dis=dis.dis;
    dis.datName=datName;
    if isfield(dis,'version')
        disOk=dis.version==version && dis.miValid;
        disp('Loading settings');
    end
end;

    dis.imageMode=1;
    dis.basePath=pwd;  %% batch mode
    dis.miIndex=0;

newDis=0;
%  disOk=0;
% if ~disOk
%     disp('Loading defaults');
%     dis=struct;
%     dis.version=version;
%     dis.datName=datName;
%     dis.basePath=pwd;
%     dis.rsMode=true;  % RSC data instead of simple image.
%     dis.readVesicleImage=false;
%     dis.readSubImage=true;
%     dis.readOnlyMode=false;
%     dis.imageMode=true;  % display micrograph; false only in StackPruner
%     dis.spMode=false;  % conventional single particle modeling
%     dis.ds=4;
%     dis.useRawAmplitudes=1;
%     dis.ccuScale=1e4;  % factor by which we scale up the raw amplitudes.
%     
%     % Box colors
%     dis.infoPath='./';
%     bColorErased=[.3 .3 .3];
%     bColorMan=[1 .7 0];  % yellow
%     bColorManNotRefined=[1 .7 0]; %orangish
%     bColorAuto=[.8 1 0];  % greenish
%     bColorBkd=[.4 1 .8]; % cyan
%     bColorVesicle=[.2 0 1]; % blue vesicle
%     bColorBadVesicle=[1 .1 .3]; % red bad vesicle
%     dis.boxColors=[bColorErased; bColorMan; bColorAuto; bColorBkd;
%         bColorVesicle; bColorManNotRefined; bColorBadVesicle];
%     dis.boxLabelColor=[0 1 0; .5 1 0; 1 1 1];
%     dis.corners=[1  1  1 1 .41 1.2 .41
%         0 .7 .7 1 .41 1.2 .41]'; %  'eroded' box corners
%     dis.lineWidth=1;
%     dis.labelFontSize=10;
%     
%     % overlay colors
%     dis.maskColor=[.9 .75 .8];
%     dis.blankColor=[1 .85 .8];
%     dis.overlapColor=[1 .83 .9];
%     dis.ghostColor=[.8 .8 1];
%     dis.ghostColorBad=[1 .7 .6];
%     dis.ghostAmp=.6;
%     
%     % initial display
%     dis.ndis=[960 960]; %%  Initial display size
%     dis.ndis=[768 768]; %%  Initial display size
%     dis.maxSize=[960 960];
%     %     dis.maxSize=dis.ndis;
%     %     dis.size=min(dis.size,dis.ndis); % Force the display to be no bigger than the image
%     dis.org=[0 0]; %pixels run from org+1 to org+size
%     dis.ax1=0;
%     dis.ax2=0;
%     dis.ax3=0;
%     
%     %     default parameters
%     dis.mode=1;  % no subtraction
%     dis.clearFigure=1;  % clear the figure when loading a new file.
%     dis.showMask=1;
%     dis.showGhosts=2;
%     dis.showBoxes=3;
%     dis.miIndex=0; % if nonzero, use this as an index into the name list.
%     dis.listParticleInfo=1;
%     dis.contrast=[4 4]; % Black, white threshold re median of normalized image
%     dis.varThresh=1000;
%     %     dis.pars=[.4 .63 1000 150 ...
%     %                 150 100 70 200 50 20];
%     dis.pars=[ 3.6   6   100  0  150  0   70  100   50   12 1 1];
%     %     pars(1:12) are:  minAmp, maxAmp; max var; rso offset;
%     %     particle blank radius, vesicle blank radius, maxBob, border, maskRadius, spect.
%     %     spectFactor, ampFactor
%     dis.pars(20)=150;  % box size in A.
%     dis.minDist=dis.pars(20)/4;  % max distance in original pixels, based on box size,
%                                  % click location from object location
%     dis.filter=[1000 20 0 0];  % inverse frequency in A; third is % inverse CTF
%     dis.infoName='';
%     dis.defaultPixA=3;
%     dis.minCCVal=1e-3;  % Value below which we figure CC is zero.
%     dis.showSpectrumInfo=false;
%     dis.zeroPreviousPicks=0;
%     dis.spectrumMaskRadiusA=100;
%     dis.spectrumScale=8;
%     dis.tFactor=1.03;  % threshold step factor
%     dis.TFactor=1.1;
%     dis.readAutopickPars=0;  % read stored autopick parameters from file
%     dis.spectrumCorrectionCoeffs=0;
%     dis.ampCorrectionCoeffs=0;
%       
%     dis.finished=0;
%     dis.miValid=0;
%     dis.goodClasses=1;
%     dis.classParticlesMode=0;
%     dis.miMode=0;
%     dis.miIndex=0;
%     newDis=1;
% end;
% 
% ok=false;
% if exist(dis.basePath,'dir')
%     cd(dis.basePath);
%     ok=true;
%     if dis.miIndex
%         disp('Loading the file miNames.mat');
%         try % check if we are pointing to an extant directory
%             ls;
%         catch
%             ok=false;
%         end;
%         if ok
%         try
%             load('miNames.mat');
%         catch
%             disp('...not found.');
%             ok=false;
%         end;
%         end;
%     end;
% end;
% if ~ok
%     dis.miIndex=0;
%     dis.miValid=0;
% end;
dis.jpegCounter=0;  % zero until a file is successfully loaded. Incremented by 'T' option.
if ~isfield(dis,'forceMicrographCoords')
    dis.forceMicrographCoords=0;
end;
%     dis.ndis=[512 512]; %%  Initial display size
dis.maxSize=[960 960];
dis.size=dis.maxSize;
dis.ndis=dis.maxSize;

dis.miValid=0;
oldB=0;

% disp('Make figure');
screenSize=get(0,'screensize');
xsiz=min(screenSize(3),dis.size(1)+3);
ysiz=min(screenSize(4)-50,dis.size(2)+3);
figure(1);
clf;
% Put the window near the top middle of the screen
set(gcf,'position',[(screenSize(3)-dis.size(1))/2 ((screenSize(4)-50)-ysiz)*0.8 xsiz ysiz],...
    'toolbar','none','menubar','none','resize','off');
% Main display
axsiz=xsiz-3;
aysiz=ysiz-3;
dis.ax1=axes('units','pixels','position',[2 3 axsiz aysiz]); %,'ticklength',[0 0]);
dis.ax2=axes('position',[.8 0 .2 .2]);
dis.ax3=axes('outerposition',[.8 0 .2 .2]);
axes(dis.ax1);
% % GetClick('init','square');
disp(['Read-only mode = ' num2str(dis.readOnlyMode)]);

picks=zeros(3,1,1,'single');
dis.classes=zeros(3,1);
ptrs=zeros(3,1,'single');

coords=zeros(1,2);
dis.currentBoxSize=dis.pars(20);
if ~isfield(dis,'autosaveJpegs'), dis.autosaveJpegs=0; end;
if ~isfield(dis,'currentFileIndex'), dis.currentFileIndex=0; end;

refreshReconstruct=0;  % flag to update the reconstruction display
previousDisMode=1;
axes3On=false;
dis.roboFitStep=0;
% dis.roboFitEndIndex=batchEnd;
% if batchEnd==0
    dis.roboFitEndIndex=inf;
% end;

% Automatic scanning, a variant of RoboFit.
scan=struct;
scan.active=false;

dis.clearFigure=1;

roboChar='naa';
b=roboChar(1);
dis.roboFitStep=1;
%%
% % interactive loop
while ((b~='q') && (b~='Q')) % q = quit; Q = quit but return to this image later
    % disp([single(b) coords]);
    %     if numel(b)>0
    switch b(1)
       case {'a' 'j' 'k' 'u' 'i'}   % autopicking
            fileOk=1;
                        b=roboChar(dis.roboFitStep);
                        dis.roboFitStep=dis.roboFitStep+1;
                    switch b(1)  % second character
                        case 'a'  % aa: Go ahead and do the auto-picking
                            disp(['Auto-picker parameters: ' num2str(dis.pars([1 2 10]))]);
                        otherwise
                            beep;
                            fileOk=0;
                     end; % 2nd character
            if fileOk
                if ~isfield(rscc,'mxCC') || max(rscc.mxCC(:))==0  % No cross correlation data
                    disp('No picking for lack of rscc preprocessor data.');
                else
                    %                     disp('Auto finding particles...');
                    % Show status in the window title bar
                    set(gcf,'name',[dis.infoName ' - autopicking...']);
                    drawnow;
                    
                    % ----- Do the autopicking here -----------------------
                    %                     Expand the mask
                    rKernel=dis.pars(9)/(dis.ds*mi.pixA);
                    netMask=(rs.mxVars>dis.pars(3))|rawMask;
                    %                     BinaryConvolve() is too slow, so use Gaussian filter
                    masks(:,:,3)=GaussFilt(netMask,.16/rKernel)>.1;
                    %                     Auto-picking
                    [coords, ovMask, rscc.mxCC2]=rspAutoParticleFinder(mi,rscc,...
                        dis,masks(:,:,3));
%                     No picking where masks(:,:,3)==1.
                    imgs(:,:,7)=150*(ovMask-.5);
                    masks(:,:,5)=ovMask>1.5;
                    imgs(:,:,8)=imscale(rscc.mxCC2,256);
                    [ptrs(3), ncf]=size(coords);
                    picks(3,1:ptrs(3),1:ncf)=coords;
                    [picks, ptrs]=rspDeleteBadAutoPicks(dis,picks,ptrs);
                    [picks, ptrs]=rspDeleteBadVesiclePicks(picks,ptrs);
                end;
                %                 update the found particles display
                refreshReconstruct=1;
                % update the mask display
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                disp([num2str(ptrs(3)) ' Particles found.']);
                set(gcf,'name',['(' num2str(dis.miIndex) ') ' dis.infoName]);
            end; % fileok

        case {'n' 'p' 'O' 'V'}  % open a file
            %             First, store the present results
                if numel(dis.infoName)>3 && dis.miValid && exist('mi','var')
                        if dis.autosaveJpegs % we save jpegs even in read-only mode.
                            CheckAndMakeDir('Picker_jpegs',1);
                            fullJpegName=['Picker_jpegs/' mi.baseFilename '_i' num2str(dis.currentFileIndex) '.jpg'];
                            print('-djpeg','-r0',fullJpegName);
                            disp(['Wrote ' fullJpegName]);
                        end;
                    if dis.readOnlyMode || dis.classParticlesMode
                        disp(['***Read-only mode*** ' num2str(sum(rspCountParticles(picks))) ' particles.']);
                    else
                        mi=rspStorePicksInMi(mi,picks,ptrs,1,dis);
                        mi.particle.autopickPars=dis.pars;
                        WriteMiFile(mi,[dis.infoName]);  % store the mi structure
                        disp([dis.infoName ' written.']);
                    end;
                    tax1=dis.ax1;
                    tax2=dis.ax2;
                    tax3=dis.ax3;
                    dis.ax1=[];
                    dis.ax2=[];
                    dis.ax3=[];
                    save(dis.datName,'dis');
                    dis.ax1=tax1;
                    dis.ax2=tax2;
                    dis.ax3=tax3;
                end;
            
            if b=='n'
                %                     disp('Get next file');
                dis.miIndex=dis.miIndex+1;
                    if dis.miIndex>numel(miNames)
                        disp('No more files!');
                        return
                    end;
                  dis.infoName=miNames{dis.miIndex};
                
            end;
            %  Actually load the files here

                [dis, mi, rscc, rs, imgs, masks, rawMask,fileOk]=rspLoadFiles2(dis,si);
                if fileOk  % successfully loaded a file
                    mi.basePath=AddSlash(pwd);
                    % Load the previous picks from the mi file.
                    %                 mi=rspUpdateMiStructure(mi);  % Change the mi.particle fields
                    if dis.zeroPreviousPicks
                        mi.particle.picks=[];
                    end;
                    [picks, ptrs, dis.classes]=rspLoadPicksFromMi(mi);
                    ptrs(4)=0; % delete old blanks
                    % % to preserve old blanks use this: [picks, ptrs, dis.classes]=rspLoadPicksFromMi(mi,picks,ptrs);
                    [picks, ptrs]=rspDeleteBadAutoPicks(dis,picks,ptrs);
                    refreshReconstruct=1;
                    partCounts=rspCountParticles(picks);
                    if any(partCounts)
                        disp([num2str(rspCountParticles(picks)) ' particles loaded']);
                    end;
                    rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                    if dis.readOnlyMode
                        set(gcf,'name',['(' num2str(dis.miIndex) ') ' mi.baseFilename ' READ-ONLY'])
                    else
                        set(gcf,'name',['(' num2str(dis.miIndex) ') ' mi.baseFilename]);
                    end;
                else
%                     disp('File couldn''t be loaded']);
                    disp('Use ''n'' to go to next, or ''O'' to start over.');
                    dis.miValid=0;
                    dis.roboFitStep=min(dis.roboFitStep,1); % skip fitting.
                end;
                       
            %         
    end;  %switch
    
    % ----------get the next click or keypress----------
    
    oldB=b(1);  % store the previous key
        dis.roboFitStep=dis.roboFitStep+1;
        if dis.roboFitStep>numel(roboChar) % wrap
            dis.roboFitStep=1;
        end;
            b=roboChar(dis.roboFitStep);
end;  % while b~='q'
%


mi=rspStorePicksInMi(mi,picks,ptrs,1,dis);
mi.particle.autopickPars=dis.pars;

if numel(dis.infoName)>3 && dis.miValid
    if dis.readOnlyMode
        disp(['Read-only mode. ' num2str(sum(rspCountParticles(picks))) ' particles.']);
    else
        WriteMiFile(mi,[dis.infoName]);  % store the mi structure
        disp([dis.infoName ' written']);
    end;
end;

dis.finished=(b=='q');  % lower-case q causes final exit.
set(gcf,'name','Done');
% if dis.finished
%     dis.miValid=0;
%     disp('Exiting');
% else
    disp('Exiting, to continue this micrograph.');
close(1);
save(datName,'dis');
if isa('si','struct')
    save([dis.basePath 'siTemp.mat'],'si');
end;
% clear batchStart

