% SimpleRSPicker.m
% F. Sigworth, Jan '13
%  See rspLoadPicksFromMi.m to see assignment of ptrs
%  See rspNewBox to see assignment of flag values
% Modified so that 'j' and 'k' have the max amp track the min amp
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

newDis=0;
%  disOk=0;
if ~disOk
    disp('Loading defaults');
    dis=struct;
    dis.version=version;
    dis.datName=datName;
    dis.basePath=pwd;
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
    dis.boxLabelColor=[0 1 0; .5 1 0; 1 1 1];
    dis.corners=[1  1  1 1 .41 1.2 .41
        0 .7 .7 1 .41 1.2 .41]'; %  'eroded' box corners
    dis.lineWidth=1;
    dis.labelFontSize=10;
    
    % overlay colors
    dis.maskColor=[.9 .75 .8];
    dis.blankColor=[1 .85 .8];
    dis.overlapColor=[1 .83 .9];
    dis.ghostColor=[.8 .8 1];
    dis.ghostColorBad=[1 .7 .6];
    dis.ghostAmp=.6;
    
    % initial display
    dis.ndis=[960 960]; %%  Initial display size
    dis.ndis=[768 768]; %%  Initial display size
    dis.maxSize=[960 960];
    %     dis.maxSize=dis.ndis;
    %     dis.size=min(dis.size,dis.ndis); % Force the display to be no bigger than the image
    dis.org=[0 0]; %pixels run from org+1 to org+size
    dis.ax1=0;
    dis.ax2=0;
    dis.ax3=0;
    
    %     default parameters
    dis.mode=1;  % no subtraction
    dis.clearFigure=1;  % clear the figure when loading a new file.
    dis.showMask=1;
    dis.showGhosts=2;
    dis.showBoxes=3;
    dis.miIndex=0; % if nonzero, use this as an index into the name list.
    dis.listParticleInfo=1;
    dis.contrast=[4 4]; % Black, white threshold re median of normalized image
    dis.varThresh=1000;
    %     dis.pars=[.4 .63 1000 150 ...
    %                 150 100 70 200 50 20];
    dis.pars=[ 3.6   6   100  0  150  0   70  100   50   12 1 1];
    %     pars(1:12) are:  minAmp, maxAmp; max var; rso offset;
    %     particle blank radius, vesicle blank radius, maxBob, border, maskRadius, spect.
    %     spectFactor, ampFactor
    dis.pars(20)=150;  % box size in A.
    dis.minDist=dis.pars(20)/4;  % max distance in original pixels, based on box size,
                                 % click location from object location
    dis.filter=[1000 20 0 0];  % inverse frequency in A; third is % inverse CTF
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
    dis.spectrumCorrectionCoeffs=0;
    dis.ampCorrectionCoeffs=0;
      
    dis.finished=0;
    dis.miValid=0;
    dis.goodClasses=1;
    dis.classParticlesMode=0;
    dis.miMode=0;
    dis.miIndex=0;
    newDis=1;
end;
%%
% dis.miIndex=1;
% dis.spectrumScale=8;
% % dis.labelFontSize=10; %%%
% % dis.zeroPreviousPicks=0;
% % dis.filter(4)=0;  % initialization
% % dis.tFactor=1.03; %%%
% % dis.pars(3)=inf; %%%
% % dis.pars(12)=1; %%%
ok=false;
if exist(dis.basePath,'dir')
    cd(dis.basePath);
    ok=true;
    if dis.miIndex
        disp('Loading the file miNames.mat');
        try
            ls;
        catch
            ok=false;
        end;
        if ok
        try
            load('miNames.mat');
        catch
            disp('...not found.');
            ok=false;
        end;
        end;
    end;
end;
if ~ok
    dis.miIndex=0;
    dis.miValid=0;
end;
dis.jpegCounter=0;  % zero until a file is successfully loaded. Incremented by 'T' option.
if ~isfield(dis,'forceMicrographCoords')
    dis.forceMicrographCoords=0;
end;
%     dis.ndis=[512 512]; %%  Initial display size
dis.maxSize=[960 960];
dis.size=dis.maxSize;
dis.ndis=dis.maxSize;
% Try to load the previous mi file

% if exist(dis.infoPath,'dir')
%     cd(dis.infoPath);
%     %     dis.miValid=0;
% end;

% Set the first command
if (dis.finished || newDis || ~dis.miValid) % We need to open a new file
    b='O';  % ask for a new file
    if dis.miIndex>1
        dis.miIndex=dis.miIndex-1; % back up to past file.
    end;
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
dis.roboFitEndIndex=inf;
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
        case {1 2 3 '.' 'w' 'W' 'x' 'v'}  % simple click or erase, info, vesicle
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
            
        case {'a' 'j' 'k' 'u' 'i'}   % autopicking
            fileOk=1;
            switch b
                case 'j' % increase amp threshold
                    dis.pars(1:2)=dis.pars(1:2)*dis.tFactor;
                case 'k' % decrease amp threshold
                    dis.pars(1:2)=dis.pars(1:2)/dis.tFactor;
                case 'i'  % increase spect limit
                    dis.pars(10)=dis.pars(10)*dis.TFactor;
                case 'u'  % decrease spect limit
                    dis.pars(10)=dis.pars(10)/dis.TFactor;
                case 'a' % two-character input for autopicking
                    if dis.roboFitStep>0
                        b=roboChar(dis.roboFitStep);
                        dis.roboFitStep=dis.roboFitStep+1;
                    else
                        [qx, qy, b]=Myginput(1,'circle'); % Put up a circle while we enter data
                        %                         [qx, qy, b]=ginput(1); % Put up a circle while we enter data
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
                            %                             dis.pars(9)=MyInput(' Mask blank radius, A',dis.pars(9));
                            
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
            end;
        case {'l' 'L'} % toggle box labels
            dis.showBoxes=dis.showBoxes+1;
            if dis.showBoxes>3
                if b=='l'
                    dis.showBoxes=2;
                else
                    dis.showBoxes=0;
                end;
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
                    %                 case 4
                    %                     disp('Vesicle display on');
            end;
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
            
            %         case {'d' 'e' 'D' 'r' 'v' }  % display mode (toggle 1-2, 2-3, 1-end,
        case {'d' 'e' 'D' 'r'}  % display mode (toggle 1-2, 2-3, 1-end,
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
                    %                 case 'v' % toggle disploay of vesicles
                    %                     if dis.mode<8
                    %                         previousDisMode=dis.mode;
                    %                         dis.mode=8;
                    %                     else
                    %                         dis.mode=previousDisMode;
                    %                     end;
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
            %             dis.filter(4)=MyInput('% Merged image sumulation',dis.filter(4));
            dis.filter(4)=0;
            
            [imgs(:,:,1:2),dis.mulr,dis.addr]=rspFilterAndScaleImages(mi,dis,rscc);
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            refreshReconstruct=1;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'g'  % toggle "ghost" vesicles
            dis.showGhosts=dis.showGhosts+1;
            if dis.showGhosts>3
                dis.showGhosts=0;
            end;
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
            %         case 'I'
            %             dis.showSpectrumInfo=~dis.showSpectrumInfo;
        case 'm'  % toggle mask display
            dis.showMask=~dis.showMask;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case {'n' 'p' 'O' 'V'}  % open a file
            %             First, store the present results
            if b~='V'
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
            end;
            
            if b=='O' % get a new file
                if ~exist('miNames','var')
                    miNames=cell(0,1);
                end;
                dis.miIndex=min(dis.miIndex,numel(miNames));
                if dis.miIndex>0
                    dis.miIndex=min(MyInput('File index ',dis.miIndex+1),numel(miNames));
                    if dis.miIndex>0
                        dis.infoName=miNames{dis.miIndex};
                    end;
                end;
                if dis.miIndex==0 % no files available.
                    disp('Find an mi file');
                    [dis.infoName,dis.infoPath] = uigetfile({'*mi.txt'},'Load Info File');
                    dis.miMode=isa(dis.infoPath,'char');
                    if dis.miMode
%                         dis.infoPath=AddSlash(dis.infoPath);
                        [dis.basePath,dis.infoPath]=ParsePath(dis.infoPath);
                        cd(dis.basePath);
                        dis.infoName=[AddSlash(dis.infoPath) dis.infoName];
%                         Make and save the mi file list
                        miNames=f2FindInfoFiles(dis.infoPath);
                        dis.miIndex=find(strcmp(dis.infoName,miNames),1);
                        dis.roboFitEndIndex=numel(miNames);
                        disp(['Saving miNames.mat']);
                        save(['miNames.mat'],'miNames');
      
                        if numel(dis.miIndex)<1 || dis.miIndex>numel(miNames)
                            disp('Setting the file index to 1.');
                            dis.miIndex=1;
                        end;
                        dis.infoName=miNames{dis.miIndex};
                    else
                        disp('Find an si file.');
                        [dis.siName,dis.basePath] = uigetfile({'*.mat'},'Load si File');
                        if ~isa(dis.siName, 'char')
                            disp('No files selected. Exiting.');
                            return
                        end;
                        si=[]; % get ready to load a new one.
                     end;
                    cd(dis.basePath);
                end;
            end; % if b~='O'
            if b=='n'
                %                     disp('Get next file');
                dis.miIndex=dis.miIndex+1;
                if dis.miIndex>numel(miNames) || ...
                        (dis.roboFitStep>0 && dis.miIndex>dis.roboFitEndIndex)
                    if dis.miIndex>numel(miNames)
                        dis.miIndex=numel(miNames); % we're at the end.
                        disp('No more files!');
                    end;
                    beep;
                    dis.roboFitStep=0;  % turn off robo fitting
                else
                    dis.infoName=miNames{dis.miIndex};
                end;
                
            elseif b=='p'
                disp('Get previous file');
                dis.miIndex=dis.miIndex-1;
                if dis.miIndex<1
                    dis.miIndex=1;
                    beep;
                    disp('Beginning file!');
                    dis.roboFitStep=0;  % turn off robo fitting
                    
                else
                    dis.infoName=miNames{dis.miIndex};
                end;
            end;
            %             else
            %
            %                 [pa, nm, ex]=fileparts(dis.infoName);
            %                 fileTypes=strcmp(ex,{'.txt';'.mat'; '.jpg'});
            %                 if any(fileTypes(1:2))
            %                     names=FindFilenames(dis.infoPath,'.+mi\....');
            %                 else
            %                     names=FindFilenames(dis.infoPath,'.+\.jpg');
            %                 end;
            %                 dis.currentFileIndex=find(strcmp(dis.infoName, names));
            %                 if numel(dis.currentFileIndex)<1
            %                     disp(['can''t find the file ' dis.infoName]);
            %
            %                     [dis.infoName,dis.infoPath] = uigetfile({'*mi.*'},'Load Info File');
            %                     if ~isa(dis.infoName,'char')
            %                         disp('No file selected');
            %                         return
            %                     end;
            %                     dis.infoPath=AddSlash(dis.infoPath);
            %                     dis.basePath=ParsePath(dis.infoPath);
            %                     cd(dis.basePath);
            %                 end;
            %
            %                 if b=='n'
            %                     %                     disp('Get next file');
            %                     ind=dis.currentFileIndex+1;
            %                     if ind>numel(names)
            %                         beep;
            %                         disp('No more files!');
            %                         dis.roboFitStep=0;  % turn off robo fitting
            %                     else
            %                         dis.infoName=names{ind};
            %                     end;
            %                     disp(' ');
            %                     disp(['File ' num2str(ind) ' of ' num2str(numel(names))]);
            %                 elseif b=='p'
            %                     disp('Get previous file');
            %                     ind=dis.currentFileIndex-1;
            %                     if ind<1
            %                         beep;
            %
            %                     else
            %                         dis.infoName=names{ind};
            %                     end;
            %                     disp(' ');
            %                     disp(['File ' num2str(ind) ' of ' num2str(numel(names))]);
            %
            %                 elseif b=='V' % reload the current file.
            %                     cd(dis.basePath);
            %                     ind=dis.currentFileIndex;
            %                 end;
            %     end;
            % end;
            %  Actually load the files here
%             if exist([dis.infoPath dis.infoName],'file') || dis.classParticlesMode  % load something
if ~dis.miMode && ~(exist('si','var') && isa(si,'struct'))
                           disp(['Loading ' dis.siName]);
                        load([dis.basePath dis.siName]);
                        nmi=numel(si.mi);
                        dis.miIndex=1;
                        dis.readOnlyMode=1;
                        dis.classParticlesMode=1;
end;



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
        case 'P'  % misc. parameters
            dis.readOnlyMode=MyInput('Read-only mode ',dis.readOnlyMode);
            dis.readAutopickPars=MyInput('Load autopick pars from file (1: no geometry; 2: all)',dis.readAutopickPars);
            dis.pars(20)=MyInput(' Box size, A ',dis.pars(20));
            dis.readVesicleImage=MyInput('Read vesicle image file ',dis.readVesicleImage);
            dis.readSubImage=MyInput('Read subtracted image file ',dis.readSubImage);
%             oldClassParticlesMode=dis.classParticlesMode;
%             dis.classParticlesMode=MyInput('Show Class Particles ',dis.classParticlesMode);
%             if dis.classParticlesMode && ~oldClassParticlesMode % Just switched it on
%                 [dis,si,miNames]=rspLoadSiClasses(dis);
%                 if ~isa(si,'struct')
%                     disp('Si class loading failed. ShowClassParticles off.');
%                     dis.classParticlesMode=0;
%                 else
%                     dis.miIndex=1

            dis.goodClasses=MyInput('Good classes ',dis.goodClasses);
%                 end;
%             end;
%             dis.miIndex=MyInput('Index into file list', dis.miIndex);
%             if dis.miIndex && ~dis.classParticlesMode
%                 if exist('Info/miNames.mat','file')
%                     load('Info/miNames.mat'); % loads miNames cell array
%                 else
%                     disp('Couldn''t find Info/miNames.mat');
%                     dis.miIndex=0;
%                 end;
%             end;
            dis.ampCorrectionCoeffs=MyInput('Amp correction coeffs ',dis.ampCorrectionCoeffs);
            dis.spectrumCorrectionCoeffs=MyInput('Spect correction coeffs ',dis.spectrumCorrectionCoeffs);
            dis.showSpectrumInfo=MyInput('Show spectrum info ',dis.showSpectrumInfo);
            dis.spectrumMaskRadiusA=MyInput('Spectrum mask radius, A ',dis.spectrumMaskRadiusA);
            %             dis.spectrumScale=MyInput('Spectrum scale-up ',dis.spectrumScale);
            dis.zeroPreviousPicks=MyInput('Zero out previous picks ',dis.zeroPreviousPicks);
            dis.autosaveJpegs=MyInput('Automatically save jpegs ', dis.autosaveJpegs);
            dis.forceMicrographCoords=MyInput('Force micrograph coords ', dis.forceMicrographCoords);
            dis.roboFitEndIndex=MyInput('Robofit end file index ',dis.roboFitEndIndex);
            %             dis.tFactor=MyInput('Amp threshold step',dis.tFactor);
            %             dis.TFactor=MyInput('Spect threshold step',dis.TFactor);
            if dis.readOnlyMode
                set(gcf,'name',[mi.baseFilename ' READ-ONLY'])
            else
                set(gcf,'name',mi.baseFilename);
            end;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'R' % Robofit, or new version of statistics collection
            % %             if b=='S' % First, run autopicking with relaxed limits
            % %                 oldPars=dis.pars;
            % %                 dis.pars(1)=0.75*dis.pars(1);
            % %                 dis.pars(10)=1.5*dis.pars(10);
            % %                 [coords, ovMask, rscc.mxCC2]=rspAutoParticleFinder(mi,rscc,...
            % %                     dis,masks(:,:,3));
            % %                 imgs(:,:,8)=imscale(rscc.mxCC2,256);
            % %                 %                     [ptrs(3), ncf]=size(coords);
            % %                 %                     picks(3,1:ptrs(3),1:ncf)=coords;
            % %                 %                 update the found particles display
            % %                 refreshReconstruct=1;
            % %                 % update the mask display
            % %                 rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            % %                 disp([num2str(ptrs(3)) ' Particles found.']);
            % %                 nParts=size(coords,1);
            % %                 %                 nParts=ptrs(3);
            % %                 stats=coords(:,[5 8]);
            % %                 %             stats=reshape(picks(3,1:nParts,[5 8]),nParts,2);
            % %                 [h1,x1]=hist(stats(:,1));
            % %                 d1=x1(2)-x1(1);
            % %                 [h2,x2]=hist(stats(:,2));
            % %                 d2=x2(2)-x2(1);
            % %                 counts=zeros(numel(x1),numel(x2));
            % %                 for i=1:numel(x1)
            % %                     for j=1:numel(x2)
            % %                         counts(i,j)=sum(stats(:,1)>=x1(i)-d1 & stats(:,2)<x2(j)+d2);
            % %                     end;
            % %                 end;
            % %                 sum(counts(:))
            % %                 q=ptrs(3)
            % %                 axes(dis.ax2);
            % %                 cla;
            % %                 axes(dis.ax3);
            % %                 contourf(x2,x1,counts);
            % %                 colorbar;
            % %                 dis.pars=oldPars;
            % %             end;
            % toggle roboFit
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
            %             if b=='S' % scan fit mode
            %                 dis.roboFitStep=2*dis.roboFitStep; % skip the first character
            % %                 If we call rspScanStep('init') with roboFitStep>0,
            % %                 this turns on scanning:
            %                 scan=rspScanStep('init',scan,dis);
            %                 disp(['Scan Fit ' str]);
            %             else
            %                 disp(['RoboFit ' str]);
            %             end;
            
            %         case 's'  % save mi file
            %             if numel(dis.infoName)>3
            %                 if dis.readOnlyMode
            %                     disp(['Read-only mode. ' num2str(sum(rspCountParticles(picks))) ' particles.']);
            %                 else
            %                     disp('Saving mi file');
            %                     mi=rspStorePicksInMi(mi,picks,ptrs,1,dis);
            %                     WriteMiFile(mi,[dis.infoPath dis.infoName]);  % store the mi structure
            %                     disp(dis.infoName);
            %                     save(dis.datName,'dis');
            %                 end;
            %             end;
            
        case 'T' % take a jpeg picture; have to click on the window (shift-click is good) to continue
            if dis.jpegCounter>0
                CheckAndMakeDir('Picker_jpegs',1);
                set(gcf,'paperpositionmode','auto');
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                fullJpegName=['Picker_jpegs/' mi.baseFilename '_' num2str(dis.jpegCounter) '.jpg'];
                print('-djpeg','-r0',fullJpegName);
                disp(['Wrote ' fullJpegName]);
                dis.jpegCounter=dis.jpegCounter+1;
                disp('Shift-click on the window to resume.');
                
            else
                disp('No data.');
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
            %             % copy the blank entries from before.
            %             ptrs(4)=oldNBlanks;
            %             picks(4,1:oldNBlanks,1:3)=oldBlanks;
            
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
            disp('x: mark a vesicle bad, removing all particles');
            disp(' ');
            disp('----Display----');
            disp('d: toggle vesicle subtraction');
            disp('l: toggle box labels');
            disp('L: cycle through box display modes');
            disp('e: show expected particle image');
            disp('r: show residual after subtraction of expected particles');
            disp('h: show histogram of particle amplitudes; <space> to clear');
            disp('g: cycle through display of "ghost" vesicles');
            disp('c: set contrast of display (low and upper histogram limits)');
            disp('f: set filtering of display');
            disp('m: toggle mask display');
            disp('w: write out particle info');
            disp('W: write out vesicle info');
            disp(' ');
            disp('----Autopicking----');
            disp('aa: run autopicking');
            disp('ap: set autopick parameters, then run');
            disp('ag: set autopick geometric parameters, then run');
            disp(['   - RSO offset: max distance inside vesicle membrane to find ' ...
                       'a particle. Zero: no limit.']);
            disp('   - Blank radius: min distance between particles');
            disp('   - Max bob: maximum distance outside a vesicle');
            disp('   - Border: min distance of particle center to edge of micrograph');
            disp('Z: erase all picks');
            disp('R: turn on robo-fitting (load next image and start autopicking)');
            disp(' ');
            disp('----Quick autopicking adjustments----');
            disp('j: raise amp threshold (more stringent)');
            disp('k: lower amp threshold (keep more particles)');
            disp('u: lower spect threshold (more stringent)');
            disp('i: raise spect threshold (keep more particles)');
            disp(' ');
            disp('----File operations----');
            disp('n: save picks and open the next file');
            disp('p: save picks and open the previous file');
            disp('O: open a new file with a file selector');
            disp('q: save picks and quit');
            disp('Q: save picks and quit, to return to this image');
            disp('T: take a picture (save a jpeg of the window in Picker_jpegs/)');
            disp('?: show this help text');
            disp(' ');
            disp('----Advanced commands----');
            disp('D: cycle through all 10 display modes');
            disp('S: robo-fitting with statistics collection, relaxed criteria');
            disp('P: set special preferences, includread-only mode');
            disp('G: set amplitude of "ghost vesicle" display');
            disp('v: mark a bad vesicle good, for re-processing');
            disp(' ');
    end;  %switch
    
    % ----------get the next click or keypress----------
    
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
%
hold off;


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
