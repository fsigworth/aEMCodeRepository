% SimpleRSPicker.m
% F. Sigworth, Jan '13
%  See rspLoadPicksFromMi.m to see assignment of ptrs
%  See rspNewBox to see assignment of flag values

version=102;

% Retrieve parameters from a file in the program directory
pa=fileparts(which('SimpleRSPicker'));
datName=[AddSlash(pa) 'SimpleRSPickerDat.mat'];
disOk=0;
if exist(datName,'file')>0
    dis=load(datName);
    dis=dis.dis;
    dis.datName=datName;
    if isfield(dis,'version')
        disOk=dis.version==version && dis.miValid;
        disp('Loading settings');
    end
end;

newDis=0;
if ~disOk
    disp('Loading defaults');
    dis=struct;
    dis.version=version;
    dis.datName=datName;
    dis.rsMode=true;  % RSC data instead of simple image.
    dis.readVesicleImage=false;
    dis.readSubImage=true;
    dis.readOnlyMode=false;
    dis.imageMode=true;  % display micrograph; false only in StackPruner
    dis.spMode=false;  % conventional single particle modeling
    dis.ds=4;
    dis.useRawAmplitudes=1;
    dis.ccuScale=1e4;  % factor by which we scale up the raw amplitudes.
    
    % Box colors
    dis.infoPath='./';
    bColorErased=[.3 .3 .3];
    bColorMan=[1 .7 0];  % yellow
    bColorManNotRefined=[1 .7 0]; %orangish
    bColorAuto=[.8 1 0];  % greenish
    bColorBkd=[.4 1 .8]; % cyan
    bColorVesicle=[.2 0 1]; % blue vesicle
    bColorBadVesicle=[1 .1 .3]; % red bad vesicle
    dis.boxColors=[bColorErased; bColorMan; bColorAuto; bColorBkd;
        bColorVesicle; bColorManNotRefined; bColorBadVesicle];
    dis.boxLabelColor=[0 1 0; .5 1 0];
    dis.corners=[1  1  1 1 .41 1.2 .41
        0 .8 .8 1 .41 1.2 .41]'; %  'eroded' box corners
    dis.lineWidth=1;
    
    % overlay colors
    dis.maskColor=[.9 .75 .8];
    dis.blankColor=[1 .85 .8];
    dis.overlapColor=[1 .83 .9];
    dis.ghostColor=[.8 .8 1];
    dis.ghostColorBad=[1 .7 .6];
    dis.ghostAmp=.6;
    
    % initial display
    dis.ndis=[960 960]; %%  Initial display size
    dis.maxSize=[960 960];
    dis.size=dis.maxSize;
    %     dis.size=min(dis.size,dis.ndis); % Force the display to be no bigger than the image
    dis.org=[0 0]; %pixels run from org+1 to org+size
    dis.ax1=0;
    dis.ax2=0;
    dis.ax3=0;
    
    %     default parameters
    dis.mode=1;  % no subtraction
    dis.clearFigure=1;  % clear the figure when loading a new file.
    dis.showMask=1;
    dis.showGhosts=1;
    dis.showBoxes=2;
    dis.miNameIndex=1; % if nonzero, use this as an index into the name list.
    dis.listParticleInfo=1;
    dis.contrast=[5 5]; % Black, white threshold re median of normalized image
    dis.varThresh=40;
%     dis.pars=[.4 .63 1000 150 ...  
%                 150 100 70 200 50 20]; 
    dis.pars=[3.6    6   1000  300  150  0   70  100   50   12 1 1];
%     pars(1:12) are:  minAmp, maxAmp; max var; rso offset;
%     particle blank radius, vesicle blank radius, maxBob, border, maskRadius, spect.
%     spectFactor, ampFactor
    dis.pars(20)=150;  % box size in A.
    dis.minDist=dis.pars(20)/5;  % distance in original pixels, based on box size
    dis.filter=[1000 20 0];  % inverse frequency in A; third is % inverse CTF
    dis.infoName='';
    dis.defaultPixA=3;
    dis.minCCVal=1e-3;  % Value below which we figure CC is zero.
    dis.showSpectrumInfo=false;
    dis.zeroPreviousPicks=0;
    dis.spectrumMaskRadiusA=100;
    dis.spectrumScale=8;
    dis.tFactor=1.03;  % threshold step factor
    dis.TFactor=1.1;
    dis.readAutopickPars=0;  % read stored autopick parameters from file
    dis.useSpectrumCorrectionTable=0;
    dis.useAmpCorrectionTable=0;
%     Tables(:,1) are defocus, (:,2) are factors to multiply spect and amp
%     values in autopicking.
    dis.spectTable=[.7  9 ; 1 10 ; 1.2 11 ; 1.5 12 ; 1.7 15 ; 1.9 16 ; 2.5 17
                    3  18 ; 3.5 20 ];
    dis.spectTable(:,2)=dis.spectTable(:,2)/16;  % use nominal value at 2um
    dis.ampTable=[1  2.0; 2 2.1; 3.5 2.6; 4 2.6; 6 2.7];
    dis.ampTable(:,2)=dis.ampTable(:,2)/2.5;

    dis.finished=0;
    dis.miValid=0;
    newDis=1;
end;
%%
% dis.miNameIndex=1;
% dis.spectrumScale=8;
dis.zeroPreviousPicks=0;
dis.filter(4)=0;  % initialization
dis.tFactor=1.03; %%%
dis.pars(3)=inf; %%%
dis.pars(12)=1; %%%
if dis.miNameIndex
    disp('Loading the file Info/allNamesSorted.mat');
    try
        load('Info/allNamesSorted.mat');
    catch
        disp('...not found.');
        dis.miNameIndex=0;
        dis.miValid=0;
    end;
end;


% Try to load the previous mi file

% if exist(dis.infoPath,'dir')
%     cd(dis.infoPath);
%     %     dis.miValid=0;
% end;

% Set the first command
if (dis.finished || newDis || ~dis.miValid) % We need to open a new file
    b='o';  % ask for a new file
else
    b='V';  % 'revert' to latest file.
end;
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

picks=zeros(1,1,3,'single');
ptrs=zeros(3,1,'single');
coords=zeros(1,2);
dis.currentBoxSize=dis.pars(20);

refreshReconstruct=0;  % flag to update the reconstruction display
previousDisMode=1;
axes3On=false;
dis.roboFitStep=0;
roboChar='naa';

% Automatic scanning, a variant of RoboFit.
scan=struct;
scan.active=false;

dis.clearFigure=1;
%%
% % interactive loop
while ((b~='q') && (b~='Q')) % q = quit; Q = quit but return to this image later
    % disp([single(b) coords]);
    %     if numel(b)>0
    switch b(1)
        case {1 2 3 '.' 'i' 'l' 'x'}  % simple click or erase, info, vesicle
            [picks, ptrs, rscc, mi, doUpdate]=rspNewBox(mi,rscc,dis,picks,ptrs,coords,b); % insert a manual pick
            refreshReconstruct=1;
            if doUpdate
                if dis.mode==8 % we're looking at vesicle models, update the image
                    %                     with the same scaling as rspLoadFiles
                    imgs(:,:,8)=dis.mulr*rscc.mVes+dis.addr;
                end;
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            end;
            if b==1 || b==2 || b=='.' || b=='x' % changed a particle
                refreshReconstruct=1;
            end;
            
        case {'a' 'j' 'k' 'K' 'J'}   % autopicking
            fileOk=1;
            switch b
                case 'j' % increase amp threshold
                    dis.pars(1)=dis.pars(1)*dis.tFactor;
                case 'k' % decrease amp threshold
                    dis.pars(1)=dis.pars(1)/dis.tFactor;
                case 'K'  % increase spect limit
                    dis.pars(10)=dis.pars(10)*dis.TFactor;
                case 'J'  % decrease spect limit
                    dis.pars(10)=dis.pars(10)/dis.TFactor;
                case 'a' % two-character input
                    if dis.roboFitStep>0
                        b=roboChar(dis.roboFitStep);
                        dis.roboFitStep=dis.roboFitStep+1;
                    else
                        [qx, qy, b]=Myginput(1,'circle'); % Put up a circle while we enter data
                    end;
                    switch b(1)  % second character
                        case 'a'  % aa: Go ahead and do the auto-picking
                            disp(['Auto-picker parameters: ' num2str(dis.pars([1 2 10]))]);
                        case 'g'  % ag: change geometry parameters, then auto-pick
                            dis.pars(4)=MyInput(' RSO offset, A  ',dis.pars(4));
                            dis.pars(5)=MyInput(' Blank radius, A',dis.pars(5));
                            dis.pars(6)=MyInput(' Overlap blanking radius, A',dis.pars(6));
                            dis.pars(7)=MyInput(' Max bob, A',dis.pars(7));
                            dis.pars(8)=MyInput(' Border, A', dis.pars(8));
                            dis.pars(9)=MyInput(' Mask blank radius, A',dis.pars(9));
                            
                        case 'p'  % ap: change parameters, then auto-pick
                            disp('Auto-picker parameters:');
                            dis.pars(1)=MyInput(' Min amplitude  ',dis.pars(1));
                            dis.pars(2)=MyInput(' Max amplitude  ',dis.pars(2));
%                             dis.pars(3)=MyInput(' Max variance   ',dis.pars(3));
                            dis.pars(10)=MyInput(' Max spectrum   ',dis.pars(10));
                            
                            % optional online parameter entry
                            %                 case 'v'  % av: set the variance threshold then go
                            %                     [dis.pars(3) ok]=rspGetGValue(dis.pars(3));
                            %                     disp(['Max variance = ' num2str(dis.pars(3))]);
                        otherwise
                            beep;
                            fileOk=0;
                    end; % 2nd character
                otherwise
                    beep;
                    fileOk=0;
            end;
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
                    [coords, ovMask, endCC]=rspAutoParticleFinder(mi,rscc,...
                        dis,masks(:,:,3));
                    imgs(:,:,7)=150*(ovMask-.5);
                    masks(:,:,5)=ovMask>1.5;
                    imgs(:,:,8)=imscale(endCC,256);
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
                set(gcf,'name',['(' num2str(dis.miNameIndex) ') ' dis.infoName]);
            end;
            
        case 'b'  % toggle box display
            dis.showBoxes=dis.showBoxes+1;
            if dis.showBoxes>4
                dis.showBoxes=0;
            end;
            dis.currentBoxSize=dis.pars(20);
            switch dis.showBoxes
                case 0
                    disp('Box display off');
                case 1
                    disp('Mini boxes');
                    dis.currentBoxSize=dis.pars(20)/2;
                case 2
                    disp('Box display on');
                case 3
                    disp('Text display on');
                case 4
                    disp('Vesicle display on');
            end;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'B'
            dis.pars(20)=MyInput(' Box size, A ',dis.pars(20));
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
        case 'c'  % set contrast
            disp('Setting contrast');
            dis.contrast(1)=MyInput('Black contrast',dis.contrast(1));
            dis.contrast(2)=MyInput('White contrast',dis.contrast(2));
            [imgs(:,:,1:2),dis.mulr,dis.addr]=rspFilterAndScaleImages(mi,dis,rscc);
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'C'  % set overlay colors
            dis.maskColor=MyInput('Mask color',dis.maskColor);
            dis.blankColor=MyInput('Blank color',dis.blankColor);
            dis.overlapColor=MyInput('Overlap color',dis.overlapColor);
            dis.ghostColor=MyInput('Ghost color',dis.ghostColor);
            dis.ghostColorBad=MyInput('Bad ghost color',dis.ghostColorBad);
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case {'d' 'e' 'D' 'r' 'v' }  % display mode (toggle 1-2, 2-3, 1-end,
            %                                 residual, vesicle fits)
            refreshReconstruct=0;
            switch b
                case 'd'
                    dis.mode=dis.mode+1;
                    if dis.mode>2-dis.spMode % lower case toggles norm / ves-subtr
                        dis.mode=1;
                    end;
                case 'e' % show modeled particles (mode 3)
                    if dis.mode~=3
                        previousDisMode=dis.mode; % remember previous mode
                        dis.mode=3;
                        refreshReconstruct=1;
                    else
                        dis.mode=previousDisMode; % toggle back to previous
                    end;
                case 'r'  % toggle display of image - particles
                    if dis.mode~=4 % We're coming from a conventional display
                        previousDisMode=dis.mode;
                        dis.mode=4;                % residual display
                        refreshReconstruct=1;
                    else
                        dis.mode=previousDisMode;  %  standard display
                    end;
                case 'v' % toggle disploay of vesicles
                    if dis.mode<8
                        previousDisMode=dis.mode;
                        dis.mode=8;
                    else
                        dis.mode=previousDisMode;
                    end;
                otherwise
                    dis.mode=dis.mode+1;
                    if dis.mode > size(imgs,3) % upper case D: cycle through all imgs
                        dis.mode=1;
                    end;
            end;
            if refreshReconstruct
                imgs(:,:,3)=rspReconstructParticles(dis,mi,picks,ptrs,rscc);
                imgs(:,:,4)=imgs(:,:,previousDisMode)-imgs(:,:,3)+dis.addr;
                refreshReconstruct=0;
            end;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'f'  % set filtering
            disp('Setting the filter:');
            dis.filter(1)=MyInput('Highpass filter, A',dis.filter(1));
            dis.filter(2)=MyInput('Lowpass filter, A',dis.filter(2));
            dis.filter(3)=MyInput('% LF compensation, percent',dis.filter(3));
            dis.filter(4)=MyInput('% Merged image sumulation',dis.filter(4));

            [imgs(:,:,1:2),dis.mulr,dis.addr]=rspFilterAndScaleImages(mi,dis,rscc);
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            refreshReconstruct=1;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'g'  % toggle "ghost" vesicles
            dis.showGhosts=~dis.showGhosts;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
        case 'G'  % set ghost amplitude
            disp('Setting ghost amplitude');
            dis.ghostAmp=MyInput('Ghost amplitude',dis.ghostAmp);
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
        case 'h'  % Toggle the autopicker histogram
            if ptrs(3)>0
                axes(dis.ax2);
                cla;
                axes(dis.ax3);
                %                 midAmp=median(mi.vesicle.s(~isnan(mi.vesicle.s)));
                amps=[];
                for i=1:ptrs(3)
                    c=squeeze(picks(3,i,:));
                    %                     if c(3)>15 && c(4)>0 % a pick
                    if c(3)>15  % a pick
                        %                             amps=[amps; c(5)*midAmp/mi.vesicle.s(c(4))];
                        amps=[amps; c(5)];
                    end;
                end;
                %                     set(ax2,'color','w');
                hist(amps,30);
                drawnow;
            end;
        case 'I'
            dis.showSpectrumInfo=~dis.showSpectrumInfo;
        case 'm'  % toggle mask display
            dis.showMask=~dis.showMask;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case {'n' 'N' 'o' 'V'}  % open a file
            %             First, store the present results
            if b~='V'
                if numel(dis.infoName)>3 && dis.miValid && exist('mi','var')
                    if dis.readOnlyMode
                        disp(['Read-only mode. ' num2str(sum(rspCountParticles(picks))) ' particles.']);
                    else
                        mi=rspStorePicksInMi(mi,picks,ptrs);
                        mi.particle.autopickPars=dis.pars;
                        WriteMiFile(mi,[dis.infoPath dis.infoName]);  % store the mi structure
                        disp([dis.infoName ' written.']);
                        disp([' with ' num2str(sum(rspCountParticles(picks))) ' particles']);
                    end;
                end;
            end;
            save(dis.datName,'dis');
            %             disp('Loading mi file');
            if b=='o' %|| ~exist([dis.infoName,dis.infoPath],'file') % get a new file
                if dis.miNameIndex && numel(allNamesSorted)>0
                    dis.miNameIndex=min(MyInput('File index ',dis.miNameIndex+1),numel(allNamesSorted));
                    if dis.miNameIndex>0
                        dis.infoName=allNamesSorted{dis.miNameIndex};
                        dis.infoPath='Info/';
                    end;
                end;
                if dis.miNameIndex<=0
                    [dis.infoName,dis.infoPath] = uigetfile({'*mi.txt'},'Load Info File');
                    if ~isa(dis.infoName,'char')
                        disp('No file selected');
                        return
                    end;
                    dis.infoPath=AddSlash(dis.infoPath);
                    dis.basePath=ParsePath(dis.infoPath);
                    cd(dis.basePath);
                end;
                
            else % open next or previous or current
                
                if dis.miNameIndex && numel(allNamesSorted)>0
                    if b=='n'
                        %                     disp('Get next file');
                        dis.miNameIndex=dis.miNameIndex+1;
                        if dis.miNameIndex>numel(allNamesSorted)
                            dis.miNameIndex=numel(allNamesSorted);
                            beep;
                            disp('No more files!');
                            dis.roboFitStep=0;  % turn off robo fitting
                        else
                            dis.infoName=allNamesSorted{dis.miNameIndex};
                        end;
                        
                    elseif b=='N'
                        disp('Get previous file');
                        dis.miNameIndex=dis.miNameIndex-1;
                        if dis.miNameIndex<1
                            dis.miNameIndex=1;
                            beep;
                            disp('Beginning file!');
                            dis.roboFitStep=0;  % turn off robo fitting
                            
                        else
                            dis.infoName=allNamesSorted{dis.miNameIndex};
                        end;
                    end;
                else
                    
                    [pa, nm, ex]=fileparts(dis.infoName);
                    fileTypes=strcmp(ex,{'.txt';'.mat'; '.jpg'});
                    if any(fileTypes(1:2))
                        names=FindFilenames(dis.infoPath,'.+mi\....');
                    else
                        names=FindFilenames(dis.infoPath,'.+\.jpg');
                    end;
                    currentIndex=find(strcmp(dis.infoName, names));
                    if numel(currentIndex)<1
                        disp(['can''t find the file ' dis.infoName]);
                        
                        [dis.infoName,dis.infoPath] = uigetfile({'*mi.*'},'Load Info File');
                        if ~isa(dis.infoName,'char')
                            disp('No file selected');
                            return
                        end;
                        dis.infoPath=AddSlash(dis.infoPath);
                        dis.basePath=ParsePath(dis.infoPath);
                        cd(dis.basePath);
                    end;
                    
                    if b=='n'
                        %                     disp('Get next file');
                        ind=currentIndex+1;
                        if ind>numel(names)
                            beep;
                            disp('No more files!');
                            dis.roboFitStep=0;  % turn off robo fitting
                        else
                            dis.infoName=names{ind};
                        end;
                        
                    elseif b=='N'
                        disp('Get previous file');
                        ind=currentIndex-1;
                        if ind<1
                            beep;
                            
                        else
                            dis.infoName=names{ind};
                        end;
                        
                    elseif b=='V' % reload the current file.
                        cd(dis.basePath);
                    end;
                end;
            end;
            %%
            if b~='z' && exist([dis.infoPath dis.infoName],'file')  % a valid file
                [dis, mi, rscc, rs, imgs, masks, rawMask,fileOk]=rspLoadFiles(dis);
                if fileOk  % successfully loaded a file
                    mi.basePath=AddSlash(pwd);
                    % Load the previous picks from the mi file.
                    %                 mi=rspUpdateMiStructure(mi);  % Change the mi.particle fields
                    if dis.zeroPreviousPicks
                        mi.particle.picks=[];
                    end;
                    [picks, ptrs]=rspLoadPicksFromMi(mi,picks,ptrs);
                    [picks, ptrs]=rspDeleteBadAutoPicks(dis,picks,ptrs);
                    refreshReconstruct=1;
                    partCounts=rspCountParticles(picks);
                    if any(partCounts)
                        disp([num2str(rspCountParticles(picks)) ' particles loaded']);
                    end;
                    rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                    if dis.readOnlyMode
                        set(gcf,'name',['(' num2str(dis.miNameIndex) ') ' mi.baseFilename ' READ-ONLY'])
                    else
                        set(gcf,'name',['(' num2str(dis.miNameIndex) ') ' mi.baseFilename]);
                    end;
                else
                    warning(['File couldn''t be loaded']);
                    dis.miValid=0;
                    dis.roboFitStep=min(dis.roboFitStep,1); % skip fitting.
                end;
            else
                dis.miValid=0;
            end;
        case 'P'  % misc. parameters
            dis.readOnlyMode=MyInput('Read-only mode ',dis.readOnlyMode);
            dis.readAutopickPars=MyInput('Load autopick pars from file ',dis.readAutopickPars);
            dis.readVesicleImage=MyInput('Read vesicle image file ',dis.readVesicleImage);
            dis.readSubImage=MyInput('Read subtracted image file ',dis.readSubImage);
            dis.miNameIndex=MyInput('Index into file list', dis.miNameIndex);
            if dis.miNameIndex
                if exist('Info/allNamesSorted.mat','file')
                    load('Info/allNamesSorted.mat');
                else
                    disp('Couldn''t find Info/allNamesSorted.mat');
                    dis.miNameIndex=0;
                end;
            end;
            dis.useAmpCorrectionTable=MyInput('Use amp correction tables ',dis.useAmpCorrectionTable);
                dis.useSpectrumCorrectionTable=dis.useAmpCorrectionTable;
            dis.showSpectrumInfo=MyInput('Show spectrum info ',dis.showSpectrumInfo);
            dis.spectrumMaskRadiusA=MyInput('Spectrum mask radius, A ',dis.spectrumMaskRadiusA);
%             dis.spectrumScale=MyInput('Spectrum scale-up ',dis.spectrumScale);
            dis.zeroPreviousPicks=MyInput('Zero out previous picks ',dis.zeroPreviousPicks);
%             dis.tFactor=MyInput('Amp threshold step',dis.tFactor);
%             dis.TFactor=MyInput('Spect threshold step',dis.TFactor);
            if dis.readOnlyMode
                set(gcf,'name',[mi.baseFilename ' READ-ONLY'])
            else
                set(gcf,'name',mi.baseFilename);
            end;
        case 'S' % new version of statistics collection
            % First, run autopicking with relaxed limits
            oldPars=dis.pars;
            dis.pars(1)=0.75*dis.pars(1);
            dis.pars(10)=1.5*dis.pars(10);            
                               [coords, ovMask, endCC]=rspAutoParticleFinder(mi,rscc,...
                        dis,masks(:,:,3));
                    imgs(:,:,8)=imscale(endCC,256);
%                     [ptrs(3), ncf]=size(coords);
%                     picks(3,1:ptrs(3),1:ncf)=coords;
                %                 update the found particles display
                refreshReconstruct=1;
                % update the mask display
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                disp([num2str(ptrs(3)) ' Particles found.']);
                nParts=size(coords,1);
%                 nParts=ptrs(3);
              stats=coords(:,[5 8]);
%             stats=reshape(picks(3,1:nParts,[5 8]),nParts,2);
            [h1,x1]=hist(stats(:,1));
            d1=x1(2)-x1(1);
            [h2,x2]=hist(stats(:,2));
            d2=x2(2)-x2(1);
            counts=zeros(numel(x1),numel(x2));
            for i=1:numel(x1)
                for j=1:numel(x2)
                    counts(i,j)=sum(stats(:,1)>=x1(i)-d1 & stats(:,2)<x2(j)+d2);
                end;
            end;
            sum(counts(:))
            q=ptrs(3)
                axes(dis.ax2);
                cla;
                axes(dis.ax3);
            contourf(x2,x1,counts);
            colorbar;
            dis.pars=oldPars;
            
        case {'R' 'S'} % toggle roboFit
            if dis.roboFitStep==0
                dis.roboFitStep=1; % turn it on.
                str='on';
                MyBusyread('init');
                dis.clearFigure=0;
            else
                dis.roboFitStep=0;
                scan.active=false; % turn of scanning too.
                str='off';
                dis.clearFigure=1;
            end;
            if b=='S' % scan fit mode
                dis.roboFitStep=2*dis.roboFitStep; % skip the first character
%                 If we call rspScanStep('init') with roboFitStep>0,
%                 this turns on scanning:
                scan=rspScanStep('init',scan,dis);
                disp(['Scan Fit ' str]);
            else
                disp(['RoboFit ' str]);
            end;

        case 's'  % save mi file
            if numel(dis.infoName)>3
                if dis.readOnlyMode
                    disp(['Read-only mode. ' num2str(sum(rspCountParticles(picks))) ' particles.']);
                else
                    disp('Saving mi file');
                    mi=rspStorePicksInMi(mi,picks,ptrs);
                    WriteMiFile(mi,[dis.infoPath dis.infoName]);  % store the mi structure
                    disp(dis.infoName);
                    save(dis.datName,'dis');
                end;
            end;
            
            %         case 'r'  % shift right
            %             disp('shift right');
            %             if dis.org(1)>=dis.ndis(1)-dis.size(1)
            %                 dis.org(1)=0;  % wrap around
            %             else
            %                 dis.org(1)=round(min(dis.org(1)+dis.size(1)/4, ...
            %                     dis.ndis(1)-dis.size(1)));
            %             end;
            %             rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
%         case 'u'  % shift up
%             disp('shift up');
%             if dis.org(2)>=dis.ndis(2)-dis.size(2)
%                 dis.org(2)=0;  % wrap around
%             else
%                 dis.org(2)=round(min(dis.org(2)+dis.size(2)/3, ...
%                     dis.ndis(2)-dis.size(2)));
%             end;
%             rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'Z'  % zero out all the boxes
            disp('Zeroing out all boxes.');
            oldNBlanks=ptrs(4);
            oldNBlanks=min(oldNBlanks,2);  %%%% force there to be max 2 blanks ------------------
            oldBlanks=picks(4,1:oldNBlanks,1:3);
            ptrs=ptrs*0;
            picks=picks*0;
            % copy the blank entries from before.
            ptrs(4)=oldNBlanks;
            picks(4,1:oldNBlanks,1:3)=oldBlanks;
            
            refreshReconstruct=1;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
        case '?'
            disp(' ');
            disp('This program is triggered by single keystrokes in the image window');
            disp(' ');
            disp('----Basic commands----');
            disp('Left button: pick a particle or erase a particle');
            disp('Right button (or ctrl-click): pick blank region');
            disp('Center button (or shift click): erase a particle');
            disp('b: change box display');
            disp('d: toggle vesicle subtraction');
            disp('n: save picks and open the next file');
            disp('N: save picks and open the previous file');
            disp('P: set special preferences, including read-only mode');
            disp('R: turn on robo-fitting (load next image and start autopicking)');
            disp('q: save picks and quit');
            disp('Z: erase all picks');
            disp('?: show this help text');
            disp(' ');
            disp('----Quick autopicking adjustments----');
            disp('j: raise amp threshold (more stringent)');
            disp('J: lower spect threshold (more stringent)');
            disp('k: lower amp threshold (keep more particles)');
            disp('K: raise spect threshold (keep more particles)');
            disp(' ');
            disp('----Advanced commands----');
            disp('D: cycle through all 10 display modes');
            disp('Q: save picks and quit, auto-restarting with this image');
            disp('e: show expected particle image');
            disp('r: show residual after subtraction of expected particles');
            disp(' ');
    end;  %switch
% -----get the next click or keypress----

    oldB=b(1);  % store the previous key
    if dis.roboFitStep==0 % not doing automatic fitting, wait for key
        [coords, b]=rspGetClick(dis);
    else
        dis.roboFitStep=dis.roboFitStep+1;
        if dis.roboFitStep>numel(roboChar) % wrap
            dis.roboFitStep=1+scan.active;
        end;
        %         RoboFit or Scan running: a keypress halts it.
        [x,y,b]=MyBusyread;
        if b>0  % a click of some sort
            dis.roboFitStep=0;  % turn off robo-fitting;
            MyBusyread('stop');
            disp('Robo-fitting interrupted.');
            [coords,b]=rspGetClick(dis);
        else  % continue robo-fitting
            b=roboChar(dis.roboFitStep);
            [scan,dis,scanDone]=rspScanStep('next',scan,dis,ptrs(3));
            if scanDone
                disp('Scan fit completed.');
                dis.roboFitStep=0;
                [coords,b]=rspGetClick(dis);
            end;
        end;
    end;
    
    if b>31
        axes(dis.ax1);  % put the main display in front.
    end;
    % end;  % if numel(b)>0
end;  % while b~='q'
%%
hold off;


mi=rspStorePicksInMi(mi,picks,ptrs);
mi.particle.autopickPars=dis.pars;

if numel(dis.infoName)>3 && dis.miValid
    if dis.readOnlyMode
        disp(['Read-only mode. ' num2str(sum(rspCountParticles(picks))) ' particles.']);
    else
        WriteMiFile(mi,[dis.infoPath dis.infoName]);  % store the mi structure
        disp([dis.infoName ' written']);
        count=rspCountParticles(picks);
        disp([' with ' num2str(count) ' particles.']);
    end;
end;

dis.finished=(b=='q');  % lower-case q causes final exit.
set(gcf,'name','Done');
if dis.finished
    %     dis.miValid=0;
    disp('Exiting');
else
    disp('Exiting, to continue this micrograph.');
end;
close(1);
save(datName,'dis');

