function varargout = Vesicle_finding_GUI(varargin)
% Find vesicles in micrographs.
% 30 Apr 20 version.
% Vesicle_finding_GUI() calls this in graphical interactive mode.
% Vesicle_finding_GUI(fileList,[contextFile]) has this run in batch mode. fileList is a
% cell array of strings being the mi.txt files to process, and the optional contextFile
% is a string giving the VFContext.mat file to use for parameters;
% otherwise the VFContextFile stored with the executable file is used.
%
% VESICLE_FINDING_GUI M-file for Vesicle_finding_GUI.fig
%      VESICLE_FINDING_GUI, by itself, creates a new VESICLE_FINDING_GUI or raises the existing
%      singleton*.
%
%      H = VESICLE_FINDING_GUI returns the handle to a new VESICLE_FINDING_GUI or the handle to
%      the existing singleton*.
%
%      VESICLE_FINDING_GUI('CALLBACK',hObject,eventData,h,...) calls the local
%      function named CALLBACK in VESICLE_FINDING_GUI.M with the given input arguments.
%
%      VESICLE_FINDING_GUI('Property','Value',...) creates a new VESICLE_FINDING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Vesicle_finding_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Vesicle_finding_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Vesicle_finding_GUI

% Last Modified by GUIDE v2.5 17-Jul-2019 15:28:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Vesicle_finding_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @Vesicle_finding_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before Vesicle_finding_GUI is made visible.
function Vesicle_finding_GUI_OpeningFcn(hObject, ~, h, varargin)

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
% varargin   command line arguments to Vesicle_finding_GUI (see VARARGIN)
%%%%% always show this main window in the center
%pixels

set( gcf, 'Units', 'pixels' );

%get your display size
% screenSize = get(0, 'ScreenSize');

%calculate the center of the display
% position = get( gcf, 'Position' );
% position(1) = (screenSize(3)-position(3))/2;
% position(2) = (screenSize(4)-position(4))/2-30;
%
% %center the window
% set( gcf, 'Position', position );

% set the button figures
set(h.pushbutton_nextName,'cdata',arrowbutton);
set(h.pushbutton_formerName,'cdata',inarrowbutton);


% initiate some data of h and some signs

% Context variables that are stored in a file
amPars.n0=256;  % Automask parameters
amPars.thresh=.1;
amPars.var=0;
amPars.width=0;
amPars.edge=.5;
amPars.denseFilt=.0015;  % A^-1
amPars.dense=.5;

sav.automaskPars=amPars;  % saved variables, stored in "VFContext.mat"
sav.fullInfoPath='';
sav.baseName='';
sav.vesicleAmps=[2e-3 5e-3 200];
sav.vesicleRadii=[100 300 10];
sav.filterFreqs=[0 .1];   % HP and LP, in A^-1
sav.contrastPars=zeros(2,1);
sav.membranePars=[1.6 60 6];
sav.beamPars=[0 0 100];  % X, Y and R
sav.black=0;
sav.white=.1;
sav.automaskFixed=0;
% sav.initTheVesicles=0;
sav.eraseOldPicks=0; % i.e don't load the old vesicles, but delete them.
sav.startPath='';
sav.fileList={};
sav.fileIndex=0;
sav.findInMask=1;

% State variables
h.sav=sav;

disp(['The current path is ' pwd]);

% Use a context file to get the parameters h.sav
% First, find our source directory
pa=fileparts(which('Vesicle_finding_GUI'));
ourContextName=[pa filesep 'VFContext.mat'];

if numel(varargin)>1 % attempt to use our given context file
    if exist(varargin{2},'file')
        ourContextName=varargin{2};
    else
        disp([varargin{2} '  (context file) not found.']);
    end;
end;

if exist(ourContextName,'file')
    disp(['Loading ' ourContextName]);
    sav1=load(ourContextName);
    %     Any field that sav1 lacks is loaded from h.sav (defaults)
    h.sav=SetDefaultValues(h.sav,sav1.sav);
    %%%%%%
    %     h.sav.initTheVesicles=0;
else
    disp(['Context file doesn''t exist:  ' ourContextName]);
end;


% Load the context file

h.doPrinting=false; %%%%%%%%
h.imageLoaded=false;
h.readOnly=false;
h.pixA=1;
h.origImage=single(0);  % Original merged image, scaled to display size
h.filtImage=single(0);  % rawImage after low and high pass filtering
h.goodVesImage=single(0);
h.badVesImage=single(0);
h.ctf=single(0);
h.filtVesImage=single(0);
h.ccVals=0;     % values of cc maxima
h.ccValsScaled=0;
% h.ccRadii=0;    % est. radius corresponding to cc maximum.
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
h.roboTrackMembranes=1;
h.makeModelVesicles=1;  % compute the model vesicles on loading
h.maxNVes=500;  % don't try to find more than this many.

% if ~isfield(h,'automaskBeamOn')
h.automaskBeamOn=0;  % don't change this if alreay defined.
% end;
pos=get(h.axes1,'position');
h.targetDisplaySize=pos(3:4);
h.displaySize=pos(3:4);
set(gcf,'WindowButtonMotionFcn',@WindowButtonMotionFcn);
set(gcf,'WindowButtonUpFcn',@WindowButtonUpFcn);

% handles for the image display
h.hSP=0;
h.hNav=0;
h.ih=0;
% h.sav.automaskPars

% Debigging output
h.debug=1;

% set the slider parameters
set(h.slider_Threshold,'Min', 0);
set(h.slider_Threshold,'Max', 1);
set(h.slider_Threshold,'SliderStep',[.1 .2]);
set(h.slider_Threshold,'Value',h.sav.automaskPars.thresh);

set(h.slider_GlobalVar,'Min', 0);
set(h.slider_GlobalVar,'Max', 1);
set(h.slider_GlobalVar,'SliderStep',[.1 .2]);
set(h.slider_GlobalVar,'Value',h.sav.automaskPars.var);

set(h.fixedCheckbox,'Value',h.sav.automaskFixed);

set(h.slider_Local,'Min', 0);
set(h.slider_Local,'Max', 1);
set(h.slider_Local,'SliderStep',[.1 .2]);
set(h.slider_Local,'Value',h.sav.automaskPars.width);

set(h.slider_Edges,'Min', 0);
set(h.slider_Edges,'Max', 1);
set(h.slider_Edges,'SliderStep',[.1 .2]);
set(h.slider_Edges,'Value',h.sav.automaskPars.edge);

set(h.slider_dense,'Min', 0);
set(h.slider_dense,'Max', 1);
set(h.slider_dense,'SliderStep',[.1 .2]);
set(h.slider_dense,'Value',h.sav.automaskPars.dense);


% set the membrane parameters
% set(h.edit_Edge,'String',sprintf('%1.1f',h.sav.membranePars(3)));
% set(h.edit_Thickness,'String',sprintf('%d',h.sav.membranePars(2)));

set(h.edit_MinRadius,'String',sprintf('%d',h.sav.vesicleRadii(1)));
set(h.edit_MaxRadius,'String',sprintf('%d',h.sav.vesicleRadii(2)));
set(h.edit_MinAmp,'String',num2str(1000*h.sav.vesicleAmps(1)));
set(h.edit_MaxAmp,'String',num2str(1000*h.sav.vesicleAmps(2)));
set(h.edit_IceAmp,'String',num2str(h.sav.vesicleAmps(3)));

set(h.edit_RadiusStep,'String',sprintf('%d',h.sav.vesicleRadii(3)));

set(h.editWhite,'String',sprintf('%4.2g',10*h.sav.white));
set(h.editBlack,'String',sprintf('%4.2g',10*h.sav.black));

set(h.edit_Highpass,'string',num2str(h.sav.filterFreqs(1)));
set(h.edit_Lowpass,'string',num2str(h.sav.filterFreqs(2)));
set(h.edit_beamX,'string',num2str(h.sav.beamPars(1)));
set(h.edit_beamY,'string',num2str(h.sav.beamPars(2)));
set(h.edit_beamR,'string',num2str(h.sav.beamPars(3)));
% h.sav.findInMask=1; %%%%
set(h.FindInMaskButton,'value',h.sav.findInMask);
% h.sav.eraseOldPicks=0; %%%%
set(h.EraseOldPicksButton,'value',h.sav.eraseOldPicks);
set(h.textFilename,'string','---');

h.output = hObject;
% h.altBasePath='';  % mi.basePath value if the original one doesn't work.
h.imageFileTypes={'m.mrc' 'mf.mrc' 'm.jpg' 'mf.jpg'};  % allowed types of
% images to load based on mi file.

set(h.MaskRadiobuttons,'SelectedObject',[]);

% Update h structure
% UIWAIT makes Vesicle_finding_GUI wait for user response (see UIRESUME)
% uiwait(h.figure1);


% % Get the figure handle, starting with the axes handle.
% fh=h.axes1;
% str=get(fh,'Type');
% while ~strcmpi(str,'figure')
%     fh=get(fh,'Parent');
%     str=get(fh,'Type');
% end;
%
% return % ###########
%
%
% % Set up the keypress function
% fhndl=@(hObject,eventdata) Vesicle_finding_GUI('KeyPressFcn',hObject,eventdata,guidata(hObject));
% set(fh,'keypressfcn',fhndl);


h.batchMode=numel(varargin)>0 && isa(varargin(1),'cell');
if h.batchMode
    disp('Batch mode.');
    h.jpegDir=AddSlash('Jpeg');
    CheckAndMakeDir(h.jpegDir,1);
    h.sav.fileList=varargin{1};
    h.sav.fileIndex=1;
else
    h=GetFileList(h);
end;
h=LoadAFile(h);
set(h.togglebutton_RoboFit,'value',h.batchMode); % start RoboFit automatically in batch
guidata(hObject, h);

if h.batchMode
    doRobofit(hObject,h);
    disp(' ');
    disp('Vesicle_finding_GUI done.');
    delete(hObject); % remove the GUI
end;
end % Vesicle_finding_GUI_OpeningFcn


% --- Outputs from this function are returned to the command line.
function varargout = Vesicle_finding_GUI_OutputFcn(~, ~, h)
%     varargout{1} = h.output;
varargout={};
end

function h=GetFileList(h);
% Check to see if we have a valid file to start with:
% In that case h.sav.startPath exists,
% [startPath h.sav.fileList{1}] is a valid filename
% and we assume that h.sav.fileIndex points to a valid file.
% Otherwise:
% We put up a file selector and pick one mi file.
% Change directory to the one above where the mi file sits,
% Also create a list of all mi
% files, which will allow moving forward and backward.
% We initialize
%   h.sav.startPath % points to the experiment directory
%   h.sav.fileList  % element is for example {'Info/thefilemi.txt'}
%   h.sav.fileIndex % points to the name you selected in the file selector.

% Check to see if we have a starting file.
%     Check if the previous file is around
if exist(h.sav.startPath,'dir') && numel(h.sav.fileList)>0 ...
        && exist([AddSlash(h.sav.startPath) h.sav.fileList{1}],'file')==2
    h.sav.fileIndex=max(1,min(h.sav.fileIndex,numel(h.sav.fileList)));
    disp(['Restarting at ' AddSlash(h.sav.startPath) h.sav.fileList{h.sav.fileIndex}]);
    cd(h.sav.startPath);
    %         [infoPath,nm,ex]=fileparts(h.sav.startFile);
    %         nm=[nm ex];
else % put up the file selector
    disp('Getting an mi.txt file');
    [nm,pa]=uigetfile('*mi.txt','Get an mi file.');
    if isnumeric(pa) % user clicked Cancel
        return
    end;
    [basePath,infoPath]=ParsePath(pa);
    h.sav.startPath=basePath;
    cd(basePath);
    % Build up the file list
    h.sav.fileList={};
    h.sav.fileIndex=1;
    d=dir(infoPath);
    for i=1:numel(d);
        if strndcmp(d(i).name,'mi.txt')
            name=[AddSlash(infoPath) d(i).name];
            h.sav.fileList=[h.sav.fileList; {name}];
            if strcmp(d(i).name,nm) % look for the one we selected
                h.sav.fileIndex=numel(h.sav.fileList);
            end;
        end;
    end;
end; % if starting file exists
end
%
%     % Check to see if we're using a star file
%     [~,~,ex]=h.sav.fileList{1};
%     h.sav.startFile=h.sav.fileList{1};
%     h.sav.fileIsStar=strcmpi(ex,'.star');
% %     If it's a star file, set h.nFiles to the number of lines in the file.
%     if h.sav.fileIsStar
%         disp(['Reading the star file: ' h.sav.fileList{1}]);
%         [h.blockNames,h.blockData]=ReadStarFile(h.sav.fileList{1},1);
%         [~,~,~,h.nFiles]=rlStarLineToMi(h.blockNames,h.blockData,0);
%     else
%

% --- Executes on button press in togglebutton_RoboFit.
function togglebutton_RoboFit_Callback(hObject, eventdata, h)
active=get(hObject,'value');
if active
    %     h.sav.fileIndex=h.sav.fileIndex+1;
    %     [h,ok]=LoadAFile(h);
    ok=true;
    if ok
        doRobofit(hObject,h);
        guidata(hObject,h);
    end;
else
    disp('----RoboFit stopping----');
end;
end

function doRobofit(hObject,h);
disp(['RoboFit starting. ' datestr(now)]);
% Robo-fit loop
%     We check the state of the button after loading each file.
while h.imageLoaded && get(h.togglebutton_RoboFit,'value')
    %         delete the old automask
    h=NewAutomask(h,false);
    set(h.togglebutton_Automask,'value',0);
    h.doTrackMembranes=0;
    %         guidata(hObject,h);
    %         Find vesicles
    %             guidata(hObject,h);
    DoFind(h.pushbutton_Find, 0, h, 0);
    h=guidata(hObject);
    pause(0.1);
    DoFind(h.pushbutton_Find, 0, h, 1);
    %     pushbutton_FindMore_Callback(h.pushbutton_FindMore,0,h);
    h=guidata(hObject);
    pause(0.1);
    %         Make a new automask
    disp('Automask on.');
    set(h.togglebutton_Automask,'value',1);
    %     guidata(hObject,h);
    togglebutton_Automask_Callback(h.togglebutton_Automask, 0, h)
    %         %
    %         %         set(h.togglebutton_Automask,'value',0);
    %         Find vesicles again
    h=guidata(hObject);
    h.doTrackMembranes=h.roboTrackMembranes;
    DoFind(h.pushbutton_Find, 0, h, 0)
    h=guidata(hObject);
    pause(0.1);
    DoFind(h.pushbutton_Find, 0, h, 1);
    
    %     pushbutton_FindMore_Callback(h.pushbutton_FindMore,0,h);
    
    %     %         Make a new automask
    %     set(h.togglebutton_Automask,'value',1);
    %     h=guidata(hObject);
    %     togglebutton_Automask_Callback(h.togglebutton_Automask, 0, h)
    %     pause(0.1);
    %     h=guidata(hObject);
    if h.batchMode && h.doPrinting
        print('-djpeg',[h.jpegDir h.mi.baseFilename 'vf.jpg']);
    end;
    %     guidata(hObject,h);
    pushbutton_nextName_Callback(h.pushbutton_nextName,0,h);
    h=guidata(hObject);
    %     h=LoadAFile(h);
end;
disp(['RoboFit ended. ' datestr(now)]);
guidata(hObject,h);
end


% --- Executes on button press in Read Only
function radiobutton25_Callback(hObject, eventdata, h)
h.readOnly=get(hObject,'Value'); % returns toggle state of radiobutton25
guidata(hObject,h);
end


function figure1_CloseRequestFcn(hObject, eventdata, h)
h=CloseFile(h);
h.miChanged=0;
rsFindVesicles3('end');  % purge the variables
pa=AddSlash(fileparts(which('Vesicle_finding_GUI')));
% Save parameters from VesiclePara.mat in the local directory
sav=h.sav;
% sav.startPath=pwd;
% sav.startFile='';
% if h.sav.fileIndex>0 && numel(h.sav.fileList)>0
%     sav.startFile=h.sav.fileList{min(h.sav.fileIndex,end)};
% end;
contextFilename='VFContext.mat';
contextArchiveName=[pa contextFilename];
disp(' ');
try
    disp(['saving ' contextFilename]);
    save(contextFilename,'sav');
catch
    disp('  ...unsuccessful!');
end;
try
    disp(['saving ' contextArchiveName]);
    save(contextArchiveName,'sav');
catch
    disp('  ...unsuccessful!');
end

disp('--Vesicle finder done.');
disp(' ');
delete(hObject);
end


function h=CloseFile(h) % store the results of operations
if ~h.imageLoaded || ~h.miChanged
    %     disp('No changes to mi file.  Not written.');
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
    h.mi.particle.picks=[];
else
    h.mi=rsMergeVesicleList2(h.mi,h.miOriginal);
end;
if numel(h.oldVesicleModel)>0  % Restore the old vesicle model
    %     h.mi.vesicleModel=h.oldVesicleModel;
end;
h.mi.log{end+1,1}=['Vesicle_finding_GUI ' TimeStamp];
outName=[h.mi.infoPath h.mi.baseFilename 'mi.txt'];
nameWritten=WriteMiFile(h.mi,outName);

disp(['Saved: ' AddSlash(pwd) nameWritten]);
disp(' ');
h.miChanged=false;
h.imageLoaded=false;

end


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% --- Executes on button press in pushbutton_LoadFile.
% function pushbutton_LoadFile_Callback(h.pushbutton_LoadFile,0,h)
% % Get a completely new file list, and load the selected file.
% [h.sav.fileList,h.sav.fileIndex]=GetFileList;
% [h,ok]=LoadFile(h);
% if ~ok
%     disp('File couldn''t be loaded.');
% end;
% guidata(hObject, h);
% end

function pushbutton_LoadFile_Callback(hObject, ~, h)
% We'll load a file list, and read the selected mi file
h=CloseFile(h);
disp('LoadFile');
h.sav.fileList={}; % Force the file selector to be used.
h=GetFileList(h);
[h,ok]=LoadAFile(h);
if ~ok
    disp('File couldn''t be loaded.');
end;
guidata(hObject, h);
end


% --- Executes on button press in pushbutton_formerName.
function pushbutton_formerName_Callback(hObject, ~, h)
h.sav.fileIndex=h.sav.fileIndex-1;
h=LoadAFile(h);
guidata(hObject,h);
end

% --- Executes on button press in pushbutton_nextName.
function pushbutton_nextName_Callback(hObject, ~, h)
h.sav.fileIndex=h.sav.fileIndex+1;
h=LoadAFile(h);
guidata(hObject,h);
end

% function [h, ok]=LoadAnotherMiFile(h,offset)
% ok=0;
% if ~h.imageLoaded
%     return
% end;
% if h.batchMode
%     h=CloseFile(h);  % save the previous one.
%     h.miChanged=0;
%     h.sav.fileIndex=h.sav.fileIndex+offset;
%     infoPath='';
%     fileName=h.sav.fileList{h.sav.fileIndex};
% else
%     infoPath=h.sav.fullInfoPath;
%     d=dir(infoPath);
%     dirIndex = [d.isdir];  % Find the index for directories
%     fileList = {d(~dirIndex).name}';  % Get a list of the files
%
%     findcurrentfile = strfind(fileList,[h.sav.baseName 'mi.']);
%     currentfile_id = find(~cellfun(@isempty,findcurrentfile));   % current working file id
%     id=currentfile_id;
%     fileFound=false;
%     while ~fileFound
%         id=id+offset;
%         if id<1 || id>numel(fileList)
%             beep;
%             return;
%         end
%         fileName=fileList{id};
%         disp(fileName);
%         fileFound=strcmp(fileName(end-3:end),'.txt');
%         h=CloseFile(h);  % save the previous one.
%         h.miChanged=0;
%     end; % while
% end;
% disp(['Reading ' fileName]);
% [mi,nameRead]=ReadMiFile([infoPath fileName]);
%
% [pa,nm,ex]=fileparts(nameRead);
% fileName=[nm ex];
% % set(h.sav.fileNameText,'String',fileName);
%
% h.sav.baseName=nm(1:numel(nm)-2);  % delete the 'mi'
% mi.basePath=AddSlash(pwd);
%
% nim=min(min(numel(mi.ctf),numel(mi.doses)),size(mi.frameSets,1));
% for i=1:nim
%     if ~isfield(mi.ctf(i),'defocus') % skip to next file, then return
%         return
%     end;
%     disp([num2str(mi.ctf(i).defocus,3) 'um.  frames: ' num2str(mi.frameSets(i,:)) '  dose: ' num2str(mi.doses(i),3)]);
% end;
%
% [h ok]=GetImageFile(mi,h); % copies mi into h.mi
% if ~ok
%     return
% end;
%
% % Update the mask
%
% % if ~isfield(mi,'mask') || ~isfield(mi.mask,'merge') || numel(mi.mask.merge)<1
% %    Old merged data: compute the base mask and insert it.
% % disp('Computing merge mask...');
% t=h.mi.mergeMatrix;
% msk=meMakeMergedImageMask(h.mi.imageSize/4,t,h.mi.imageSize/(4*h.borderFraction));
% h.mi=meInsertMask(msk,h.mi,1);
% % end;
% %     Point to the end of the stack
% h.maskIndex=numel(h.mi.mask);
%
%
% h=InitDisplay(h);
%
% end
%



function [h,ok]=LoadAFile(h)  % load the mi file with name h.sav.fileList{h.sav.fileIndex}
h=CloseFile(h);
ok=false;
h.imageLoaded=false;
if h.sav.fileIndex<1 || h.sav.fileIndex>numel(h.sav.fileList)
    beep;
    disp('No file loaded.');
    return
else
    fileName=h.sav.fileList{h.sav.fileIndex};
    if ~exist(fileName,'file');
        disp(['File not found: ' fileName]);
        return
    end;
end;
disp(['Reading ' fileName]);
mi=ReadMiFile(fileName);  % loads mi
mi.infoPath=AddSlash(fileparts(fileName)); % extract the path for saving later.
nim=max(1,min(min(numel(mi.ctf),numel(mi.doses)),size(mi.frameSets,1)));
for i=1:nim
    disp([num2str(mi.ctf(i).defocus,3) 'um.  dose: ' num2str(mi.doses(i),3)]);
end;
h.sav.baseName=mi.baseFilename;
mi.ctf.B=100;


[h, ok]=GetImageFile(mi,h);
if ok
    h.imageLoaded=true;
    h=InitDisplay(h);
else
    disp('Image couldn''t be loaded');
end;
end

% function mis=ZeroOutParticles(mis)
% mis.particle=struct;
% mis.particle.picks=[];
% mis.particle.autopickPars=[];
% end
%
%
%
% function mis=ZeroOutVesicles(mis,h)
% % Initialize all the vesicle fields in the mi structure
% % membrane model
% vLipid=h.sav.membranePars(1);
% thk=h.sav.membranePars(2);
% rise=h.sav.membranePars(3);
% pixA=mis.pixA;
% % Create the model, which is sampled in units of the original pixel size.
% nm0=ceil(thk/(2*pixA))*2+3;  % array for vesicle model; 60A nominal
% mis.vesicleModel=fuzzymask(nm0,1,thk/pixA/2,rise/pixA)...
%     *vLipid;  % units of V
% % Zero out the previous picks
% miTemp=meCreateMicrographInfoStruct14;
% mis.vesicle=miTemp.vesicle;
% % mis.vesicle.x=[];
% % mis.vesicle.y=[];
% % mis.vesicle.r=[];
% % mis.vesicle.s=[];
% % mis.vesicle.ok=[];
% % mis.vesicle.shiftX=[];
% % mis.vesicle.shiftY=[];
% % mis.vesicle.shiftOk=[];
% % mis.vesicle.af=[];
% % mis.vesicle.refined=0;
% % mis.vesicle.extraPeaks=[];
% % mis.vesicle.extraSD=0;
% % mis.vesicle.extraS=[];
% end

function [h,ok]=GetImageFile(mi,h)
% Having loaded an mi file, load the image and initialize variables.

set(h.textFilename,'string','---');
% set(h.togglebutton_Beam,'value',0);
set(h.togglebutton_Automask,'value',0);
set(h.togglebutton_PaintMask,'value',0);
h.imageLoaded=false;

% update the vesicle info to the canonical nv x 4 logical matrix
nv=numel(mi.vesicle.x);
if ~(isfield(mi.vesicle,'ok')...
        && size(mi.vesicle.ok,1)==nv)
    mi.vesicle.ok=true(nv,1);
end;
if size(mi.vesicle.ok,2)<4 && nv>0
    mi.vesicle.ok(nv,4)=false;  % extend the array
end;

% Now set h.displaySize to the nearest nice numbers to the target display size
sz=zeros(2,2);
sz(1,:)=NextNiceNumber(h.targetDisplaySize,5,4)-h.targetDisplaySize;
sz(2,:)=NextNiceNumber(h.targetDisplaySize,5,-4)-h.targetDisplaySize;
[~,best]=min(sum(sz.^2,2));
h.displaySize=sz(best,:)+h.targetDisplaySize;


% Load the merged image
sz=max(h.displaySize);
% targetSizes=[2/3*sz sz];
[h.origImage,M,ok]=meLoadNormalizedImage(mi,h.displaySize,'m');
if ok
    % imageBasename=[mi.procPath mi.baseFilename 'm.mrc'];
    % sufExts={'s.mrc' 'z.tif' '.mrc'};
    % [fullMergedImageName,ok]=CheckForAltImage(imageBasename,sufExts);  % valid filename?  Load it
    %
    % if ok % Read the merged image
    %     m=single(ReadEMFile(fullMergedImageName));
    %     disp([fullMergedImageName ' loaded.']);
    % else
    %     % try for reading the raw micrograph. We then subtract the median and
    %     %     scale it to reflect fractional image intensity
    %     fullImageName=[mi.imagePath mi.imageFilenames{1}];
    %     if exist(fullImageName,'file') ...
    %             && isfield(mi,'imageNormScale') && mi.imageNormScale~=0
    %         m=single(ReadEMFile(fullImageName));
    %         if ~isfield(mi,'imageMedian')
    %             mi.imageMedian=median(m(:));
    %         end;
    %         m=(m-mi.imageMedian)*mi.imageNormScale;
    %         ok=true;
    %         disp([fullImageName ' (micrograph) loaded.']);
    %         disp(['  imageNormScale = ' num2str(mi.imageNormScale)]);
    %     elseif ~isfield(mi,'imageNormScale') || mi.imageNormScale==0
    %         disp([fullImageName ' (micrograph) was found...']);
    %         disp('  but mi.imageNormScale isn''t present.');
    %     end;
    % end;
    % if ok
    %     %     Determine what downsampling might already have occurred
    %     szm=size(m);
    %     % estimate the downsampling factor and origin shift from the micrograph to
    %     % our image.
    %     if any(mi.imageSize==0) % No image size given, we'll assume this is the
    % %         original micrograph size.
    %         mi.imageSize=size(m);
    %         disp(['mi.imageSize assigned to be ' num2str(mi.imageSize)]);
    %     end;
    %
    %     dsMicrograph=round(max(mi.imageSize./szm)); % we assume that any cropping
    %     %                                               or padding is small.
    %     h.M0=meGetImageScaling(mi.imageSize,szm,dsMicrograph);
    %     %
    
    %         Get our working "origImage" of size h.displaySize, and the
    %         affine matrix from its coordinates to the micrograph
    %         coordinates.
    %     [M1,h.origImage]=meGetImageScaling(m,h.displaySize);
    h.M2=M; % Overall image scaling matrix
    h.ds0=h.M2(1,1); % Overall downsampling factor
    %     h.ds0Shift=-h.M2(1:2,3)'; % Shift in micrograph coordinates
    %     micrographXY = M2*localXY (assuming zero-origin arrays) or
    %     micrographXY = ds0*localXY-ds0Shift;
    h.pixA=h.ds0*mi.pixA; % scale of displayed image.
    ok=true;
    
    h.oldFilterFreqs=[0 0];
    set(h.textFilename,'string',mi.baseFilename);
    h.e1ImageOffset=0; % changed, was =5.
    e1Image=h.origImage+h.e1ImageOffset;  % offset it from zero so we
    %         can analyze its distribution.
    h.exp1Image=Downsample(e1Image,h.displaySize/2); % displaySize is a multiple of 4.
    %     We'll also make masks of h.displaySize/2.
    h.imageLoaded=true;
    h.rawImage=single(0);
    h.automaskOn=false;
    set(h.togglebutton_Automask,'value',false);
    %     Clear all the processed images.
    h.filtImage=single(0);
    h.rawVesImage=single(0);
    h.filtVesImage=single(0);
    %     cla(h.axes3);  % trying to get around bug in Matlab 2015a
    
    % initialize the mask
    t=mi.mergeMatrix;
    baseMaskN=round(mi.imageSize/(2*h.ds0));
    msk=meMakeMergedImageMask(baseMaskN,t,baseMaskN/h.borderFraction);
    mi=InsertOurMask(msk,mi,h.ds0*2,1,'AND');
    %     Remember, when we retrieve our mask we must pad it to h.displaySize
    %     This is done by GetOurMask().
    
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
    %     h.ccRadii=0;
    h.markedVesicleIndex=0;
    h.ifImageComp=zeros(256,'single');
    %     h.ifImageFlat=zeros(256,'single');
    %     if h.maskIndex<3 % no audomasking has been done
    %         set(h.togglebutton_Automask,'value',true);
    %     end;
else
    disp('Couldn''t find the image file.');
end
end

% ----InitDisplay----
function h=InitDisplay(h)
if ~h.imageLoaded
    return
end;
% set(h.figure1,'pointer','watch');
% osize=size(h.origImage); % loaded image


%         Create the displayed image
% h.rawImage=DownsampleGeneral(h.origImage,h.displaySize);

h=CreateE1Map(h); %Create the automask image.

if h.maskIndex<3  % Turn off automask if it hasn't been done before.
    set(h.togglebutton_Automask,'value',false);
end;

if h.makeModelVesicles
    %     Handle old mi versions :(
    if ~isfield(h.mi.vesicle,'ok') || size(h.mi.vesicle.ok,1)<1
        h.mi.vesicle.ok=false;  % make one entry.
    else % something these
        nv=size(h.mi.vesicle.ok,1);
        disp(['Making ' num2str(nv) ' model vesicle' ess(nv)]);
    end;
    ne=size(h.mi.vesicle.ok,2);
    if ne<4
        h.mi.vesicle.ok(1,4)=false;  % extend it.
    end;
    
    goodVes=all(h.mi.vesicle.ok(:,1:2),2); % vesicles in range
    badVes=(h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2)); % found, but not in range
    scl.n=h.displaySize;
    scl.M=h.M2;
    h.goodVesImage=meMakeModelVesicles(h.mi,scl,find(goodVes));
    h.badVesImage=meMakeModelVesicles(h.mi,scl,find(badVes));
    %     disp('...done');
end;
h.ctf=meGetEffectiveCTF(h.mi,h.displaySize,h.ds0);
h=UpdateDisplayFiltering(h);
h=UpdateAutomaskBeam(h);  % Also calls ShowImage.


% Get the figure handle, starting with the axes handle.
fh=h.axes1;
str=get(fh,'Type');
while ~strcmpi(str,'figure')
    fh=get(fh,'Parent');
    str=get(fh,'Type');
end;

% Set up the keypress function
fhndl=@(hObject,eventdata) Vesicle_finding_GUI('KeyPressFcn',hObject,eventdata,guidata(hObject));
set(fh,'keypressfcn',fhndl);
end


% UpdateDisplayFiltering
function h=UpdateDisplayFiltering(h,doRawImage)
if nargin<2
    doRawImage=1;
end;
h.filtImage=h.origImage;
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
    
    % Draw the scatterplot of existing vesicles **redundant with code in
    % doFind()
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
            xlabel('Radius, A');
            ylabel('Amplitude');
            title(h.axes3,[num2str(sum(goodVes)) ' vesicles in range']);
        else
            title(h.axes3,' '); % clear the title
            goodVes=false;
            badVes=false;
        end;
    end;
    nv=numel(h.mi.vesicle.x);
    
    % if numel(h.filtImage)>0  % an image has been loaded
    %     msk=rot90(meGetMask(h.mi,size(h.filtImage),1:h.maskIndex));
    
    % Get a mask at full display size.
    msk=rot90(GetOurMask(h.mi,size(h.filtImage),h.ds0,1:h.maskIndex));
    % Set the display mode
    showMask=1;
    showAmps=1;
    showCircles=1;
    %     showGhosts=0;
    imData=h.filtImage;
    switch h.displayMode
        case 0 % default: image and circles.
        case 1 % subtracted with circles
            imData=h.filtImage-h.filtVesImage;
        case 2 % subtracted only
            imData=h.filtImage-h.filtVesImage;
            %             showGhosts=h.makeModelVesicles;
            showCircles=0;
            showAmps=0;
        case 3 % Unsubtracted
            showCircles=0;
            showAmps=0;
            %                 imData=h.ccValsScaled(1:n(1),1:n(2));
            %                 showGhosts=h.makeModelVesicles;
            %             else
            %                 return
            %             end;
            %         case 3
            % %             Show CTF-corrected image.
            %             imData=Downsample(h.ifImageComp,size(h.filtImage));
            % %             showGhosts=h.makeModelVesicles;
            %             showCircles=0;
    end;
    %     theImage =  repmat(rot90(imscale(imData,256,1e-3)),[1 1 3]);
    %     midValue=(h.e1CtrValue-h.e1ImageOffset)/(h.mi.doses(1)*h.mi.cpe)
    midValue=0;
    %%
    theImage =  repmat(rot90(256*(imData-midValue-h.sav.black)/(h.sav.white-h.sav.black)+128),[1 1 3]);
    nx=size(h.origImage,1);
    ny=size(h.origImage,2);
    
    %     if showGhosts
    %         ghostColor=[.7 .7 1];
    %         ghostColorBad=[1 .5 .35];
    %         dotW=2;
    %
    %         xCtrs=max(dotW+1,min(nx-dotW,round(h.mi.vesicle.x/h.ds0+1)));
    %         yCtrs=max(dotW+1,min(ny-dotW,round(h.mi.vesicle.y/h.ds0+1)));
    %         if any(goodVes)
    %             gves=imscale(max(-h.goodVesImage,-1e-9),1);
    %             for i=find(goodVes)'
    %                 gves(xCtrs(i)-dotW:xCtrs(i)+dotW,yCtrs(i)-dotW:yCtrs(i)+dotW)=1.5;
    %             end;
    %         else
    %             gves=0;
    %         end;
    %         if any(badVes)
    %             bves=imscale(max(-h.badVesImage,-1e-9),1);
    %             for i=find(badVes)'
    %                 bves(xCtrs(i)-dotW:xCtrs(i)+dotW,yCtrs(i)-dotW:yCtrs(i)+dotW)=1.5;
    %             end;
    %
    %             %             hp1=plot(h.mi.vesicle.x(goodVes)/h.ds0+1,...
    %             %             ny-(h.mi.vesicle.y(goodVes)/h.ds0),'b.','markersize',10);
    %
    %         else
    %             bves=0;
    %         end;
    %
    %         color=(1-ghostColor); % color to subtract for membrane
    %         for i=1:3
    %             theImage(:,:,i)=theImage(:,:,i).*(1-rot90(gves)*color(i));
    %         end;
    %         color=(1-ghostColorBad); % color for bad vesicle
    %         for i=1:3
    %             theImage(:,:,i)=theImage(:,:,i).*(1-rot90(bves)*color(i));
    %         end;
    %     end;
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
    %         ny=size(h.origImage,2);
    %         hp1=plot(h.mi.vesicle.x(goodVes)/h.ds0+1,...
    %             ny-(h.mi.vesicle.y(goodVes)/h.ds0),'b.','markersize',10);
    %         hp2=plot(h.mi.vesicle.x(badVes)/h.ds0+1,...
    %             ny-(h.mi.vesicle.y(badVes)/h.ds0),'r.','markersize',10);
    %         hold off;
    %         set(hp1,'HitTest','off');
    %         set(hp2,'HitTest','off');
    %     end;
    nv=numel(h.mi.vesicle.x);
    if nv>0 % something to show
        globalXY=[h.mi.vesicle.x h.mi.vesicle.y ones(nv,1,'single')]';
        vesXY=double(1+h.M2\globalXY); % convert to window coordinates
        if showAmps
            for i=1:nv
                x=vesXY(1,i);
                y=ny-vesXY(2,i);
                %             double((h.mi.vesicle.x(i)-h.M2(1,3))/h.ds0+1);
                %             y=double(ny-(h.mi.vesicle.y(i)-h.M2(2,3))/h.ds0+1);
                amp=h.mi.vesicle.s(i,1)*1000;
                text(x,y,num2str(amp,2),'color',[.8 .8 0],'fontsize',12,'PickableParts','none');
            end;
        end;
        
        if showCircles
            hold on
            for i=1:nv
                r1=h.mi.vesicle.r(i,:)/h.ds0;
                if r1(1)<1
                    continue
                end;
                [x,y]=CircleLineSegments(r1,min(10,100/r1(1)));
                x=x+vesXY(1,i);
                y=y+ny-vesXY(2,i);
                if goodVes(i)
                    plot(x,y,'b-','HitTest','off');
                elseif badVes(i)
                    plot(x,y,'r-','HitTest','off');
                end;
            end;
            hold off
        end;
        
    end;
    
    fhndl=@(hObject,eventdata) Vesicle_finding_GUI('axes1_ButtonDownFcn',hObject,eventdata,guidata(hObject));
    set(ih,'buttondownfcn',fhndl);
    set(ih,'HitTest','on');
    h.ih=ih;
    h.theImage=uint8(theImage);
    
end
end

% Vesicle_finding_GUI('axes1_ButtonDownFcn',hObject,eventdata,guidata(hObject))


function editWhite_Callback(hObject, eventdata, h)
q=str2double(get(hObject,'String'));
if ~isnan(q)
    h.sav.white=q/10;
    set(hObject,'string',sprintf('%4.2g',q));
    h=ShowImage(h);
    guidata(hObject, h);
end;
end


function editBlack_Callback(hObject, eventdata, h)
q=str2double(get(hObject,'String'));
if ~isnan(q)
    h.sav.black=q/10;
    set(hObject,'string',sprintf('%4.2g',q));
    h=ShowImage(h);
    guidata(hObject, h);
end;
end


% --- Executes on button press in AutoContrast.
function pushbutton_AutoContrast_Callback(hObject, eventdata, h)
% theImage =  repmat(rot90(256*(imData-midValue-h.sav.black)/(h.sav.white-h.sav.black)+128),[1 1 3]);
[~,mulr,addr]=imscale(h.filtImage,1,[4.5 5.5]);
h.sav.black=(.5-addr)/mulr;
h.sav.white=(1.5-addr)/mulr;
set(h.editWhite,'String',sprintf('%4.2g',10*h.sav.white));
set(h.editBlack,'String',sprintf('%4.2g',10*h.sav.black));
h=ShowImage(h);
guidata(hObject,h);
end



%%%%%%%%%%%%%%%%%%%%  Filter frequencies

function edit_Highpass_Callback(hObject, eventdata, h)
% Hints: get(hObject,'String') returns contents of edit_Highpass as text
q=str2double(get(hObject,'String'));
if ~isnan(q)
    if q<1e-6
        q=0;
    end;
    h.sav.filterFreqs(1)=q;
    set(hObject,'string',num2str(q));
    h=UpdateDisplayFiltering(h);
    h=ShowImage(h);
    guidata(hObject, h);
end;
end


function edit_Lowpass_Callback(hObject, eventdata, h)
q=str2double(get(hObject,'String'));
if ~isnan(q)
    if q<1e-4 || q>1
        q=0;
    end;
    h.sav.filterFreqs(2)=q;
    set(hObject,'string',num2str(q));
    h=UpdateDisplayFiltering(h);
    h=ShowImage(h);
    guidata(hObject, h);
    
end
end

% --- Executes on button press in pushbutton_AdjustContrast.
function pushbutton_AdjustContrast_Callback(~, ~, h)
if h.ih  % The image has been initialized
    imcontrast(h.axes1);
end
end


%%%%%%%%%%%%% MANUAL MASK %%%%%%%%%%%%%
% Local functions to handle padded, downsampled display images.

function mi=InsertOurMask(msk,mi,ds0,index,mode)
n=round(mi.imageSize/ds0); % actual size that will be inserted
newMask=Crop(msk,n); % we assume this part is centered in msk
mi=meInsertMask(newMask,mi,index,mode);
end

function msk=GetOurMask(mi,n,ds0,indices)
% Get a mask at downsample ratio ds0 (from mi.imageSize) and padded to n.
n0=round(mi.padImageSize/ds0); % expected size of the stored mask
localMask=meGetMask(mi,n0,indices);
msk=logical(Crop(localMask,n));
end

% --- Executes on button press in pushbutton_OutlineMask.
function pushbutton_OutlineMask_Callback(hObject, eventdata, h)
% Use the stored property 'value' in the axes1 button down function to
% actually initiate the manual mask input.
active=get(hObject,'value');
h.outlineMaskActive=active;
if ~active
    return
end;
if ~h.imageLoaded
    set(h.pushbutton_OutlineMask,'value',0);
elseif ~isfield(h.mi,'mask')
    h.mi.mask=struct('merge',[],'encoding',[],'data',[]);
end;
q=get(h.pushbutton_OutlineMask,'value');
if h.imageLoaded && q % We're doing manual masking
    %         try
    hFH = imfreehand();
    %     binaryImage =~rot90(hFH.createMask(h.ih),3);
    binaryImage=createMask(hFH,h.ih);
    delete(hFH);  % needed for bug in v2015a prerelease
    binaryImage =~rot90(binaryImage,3);
    h.maskIndex=max(3,h.maskIndex)+1;
    %     h.mi=meInsertMask(binaryImage,h.mi,h.maskIndex,'AND');
    h.mi=InsertOurMask(binaryImage,h.mi,h.ds0,h.maskIndex,'AND');
    
    h.miChanged=1;
    h=ShowImage(h);
    %         catch   % Error occurred, exit manual mask mode.
    %  disp('error in freehand');
    %             set(h.pushbutton_OutlineMask,'value',0);
    %             return
    %         end;
    
    % Create a binary image ("mask") from the ROI object.
    set(h.pushbutton_OutlineMask,'value',0);
    h.outlineMaskActive=0;
end;
guidata(hObject, h);
end

% --- Executes on button press in togglebutton_PaintMask.
function togglebutton_PaintMask_Callback(hObject, eventdata, h)
active=get(hObject,'value');
h.manualMaskActive=active;
guidata(hObject,h);
end



function h=InsertPaintedMask(h)
if numel(h.manualMaskCoords)<2  % nothing drawn, remain active.
    return
end;
% We're done painting a mask.  First turn off the button
set(h.togglebutton_PaintMask,'value',0);
h.manualMaskActive=0;
n=round(size(h.origImage)/2);  % make a half-sized mask
msk=false(n);
r=h.manualMaskDiameter/4;
coords=round((h.manualMaskCoords+1)/2);
h.manualMaskCoords=zeros(0,2,'single');
npts=size(coords,1);

ru=ceil(r*2+1);
template=single(disc(ru,r));

for i=1:npts
    msk1=ExtractImage(template,coords(i,:),n,1);
    msk=msk | msk1;
end;
msk=fliplr(msk);
h.maskIndex=max(3,h.maskIndex)+1;
% h.mi=meInsertMask(~msk,h.mi,h.maskIndex,'AND');
h.mi=InsertOurMask(~msk,h.mi,h.ds0*2,h.maskIndex,'AND');
h.miChanged=1;
h=ShowImage(h);
end


% --- Executes on button press in pushbutton_Invert.
function pushbutton_Invert_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton_Invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
msk=meGetMask(h.mi,0,h.maskIndex);  % get mask at native size.
h.mi=meInsertMask(~msk,h.mi,h.maskIndex,...
    h.mi.mask(h.maskIndex).merge);
h.miChanged=1;
h=ShowImage(h);
guidata(hObject, h);
end


% --- Executes on button press in pushbutton_Intersect.
function pushbutton_Intersect_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton_Intersect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
switch h.mi.mask(h.maskIndex).merge
    case 'AND'
        h.mi.mask(h.maskIndex).merge='OR';
    case 'OR'
        h.mi.mask(h.maskIndex).merge='AND';
end;
h.miChanged=1;
h=ShowImage(h);
guidata(hObject, h);
end

% --- Executes on button press in pushbutton_Undo.
function pushbutton_Undo_Callback(hObject, eventdata, h)
h.maskIndex=max(1,h.maskIndex-1);
h=ShowImage(h);
guidata(hObject, h);
end


% --- Executes on button press in pushbutton_Redo.
function pushbutton_Redo_Callback(hObject, eventdata, h)
h.maskIndex=min(h.maskIndex+1,numel(h.mi.mask));
h=ShowImage(h);
guidata(hObject, h);
end


% ---- Executes on change of manual mask radiobuttons
function MaskRadiobuttonsCallback(hObject, eventdata, h)
p=get(hObject,'SelectedObject');
val=get(p,'string');
h=guidata(hObject);  % We'll pick up first 3 digits of the string.
h.manualMaskDiameter=str2double(val(1:3))/h.pixA;
guidata(hObject,h);
end

% --- Executes during object creation, after setting all properties.
function MaskRadiobuttons_CreateFcn(hObject, eventdata, h)
set(hObject,'selectionchangefcn',@MaskRadiobuttonsCallback);
set(hObject,'SelectedObject',[]);
%         h.manualMaskDiameter=0;
%
%         children=get(hObject,'Children');
%         set(hObject,'SelectedObject',children(1));
%         set(children(1),'value',1);
%         q=get(children(1),'string');
%         q=q(1:3)
% %         h.manualMaskDiameter=str2double(q)/h.pixA;
%         guidata(hObject,h);
end




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



function h=CreateE1Map(h)  % Create the e1Map, based on the 1st exposure
%                           (if available) for automasking.
%   Create the global "dense" mask
if numel(h.exp1Image)>1 && ~any(isnan(h.exp1Image(:))) % Something there
    %     filter it, make a histogram, set the threshold.  We take the
    %     reference intensity of the "clear" part of the micrograph to be
    %     the upper half-maximum point of the histogram.
    M2Inv=inv(h.M2);
    
    exp1Mask=ones(round(h.mi.padImageSize*M2Inv(1,1)/2),'single'); % useful part of the image
    exp1Mask=logical(Crop(exp1Mask,size(h.exp1Image))); % expand to the image's size
    exp1Mask=GetOurMask(h.mi,size(h.exp1Image),h.M2(1,1)*2,1);
    exp1f=GaussFiltDCT(h.exp1Image,h.sav.automaskPars.denseFilt*h.pixA*2);
    exp1fSel=exp1f(exp1Mask);
    med=abs(median(exp1fSel(:)));  % make sure it's positive, should be.
%     bins=.75*med:med/200:1.25*med;
    [e1Hist,bins]=hist(exp1fSel(:),50);
    [~,e1Peak]=max(e1Hist);
    e1Upper=find(e1Hist>e1Peak/3,1,'last');
    e1Lower=find(e1Hist>e1Peak/3,1);
    %     sigma5=5*(e1Upper-e1Mode+.5);
    %     e1Ctr=bins(e1Mode);
    e1Mode=bins(e1Peak);
%     e1Ctr=bins(e1Upper);
    e1Ctr=bins(e1Lower);
    if h.sav.automaskFixed
        e1Ctr=h.e1CtrValue;
    else
        h.e1CtrValue=e1Ctr;
    end;
    
    %     sigma1=std(exp1f(abs(exp1f-e1Ctr)<sigma5));
    %     e1mx=e1Ctr+sigma1;
    %     e1mn=e1Ctr-10*sigma1;
    %     h.e1Map=(exp1f-e1mn)/(e1mx-e1mn);
%     h.e1Map=(10*(exp1f/e1Ctr-.9)); % expand the top 10%
    h.e1Map=.5*((exp1f-e1Ctr)/(e1Mode-e1Ctr)+1); % expand upper half to 0..1
else
    h.e1Map=ones(size(h.exp1Image),'single');  % all ones
end;

% Also create the varMap.
f1=(h.pixA*2)/20;  % lowpass cutoff of bp filter for variance
f2=(h.pixA*2)/60;  % highpass cutoff
varFilt=(h.pixA*2)/400; % filtering of the variance
workImage=DownsampleGeneral(h.origImage,h.displaySize/2);
varMap=GaussFiltDCT((SharpHP(SharpFilt(workImage,f1),f2)).^2,varFilt);
varMapPts=varMap(exp1Mask);
mxVal=max(varMapPts(:));  % make sure it's positive, should be.
% bins=0:mxVal/100:mxVal;
% varHist=hist(varMap(:),bins);
% [~,modeIndex]=max(varHist);
% varMode=bins(modeIndex);
% h.varMap=varMap/(10*varMode);
h.varMap=varMap/mxVal;
end

function h=InitAutomask(h)
h.maskIndex=max(h.maskIndex,3);  % Show this mask
m=h.origImage-h.goodVesImage-h.badVesImage;
m=Downsample(m,size(m)/2);
h.ifImage=meCTFInverseFilter(m,h.mi,1,0,0);  % totally inverse filtered
% h.ifImageFlat=GaussFilt(meCTFInverseFilter(m,h.mi,1,0,.0005),.05);
h.ifImageComp=GaussFilt(meCTFInverseFilter(m,h.mi,1,.0005,.0005),.1);

% q=gcf;
% figure(2);
% subplot(221); imags(h.ifImageFlat);
% subplot(222); imags(h.ifImageComp);
% figure(q);
% pwd
% save('h.mat','h');

%
% h.automaskLims=zeros(1,4);
% % h.automaskLims(1:2)=Place(h.ifImageFlat(:),[.01 .8]);
% h.automaskLims(3:4)=Place(h.ifImageComp(:),[.01 .8]);
if ~isfield(h.mi,'mask')
    h.mi.mask=struct('merge',[],'encoding',[],'data',[]);
end;
end

function vals=Place(m,pctls)
np=numel(m);
[sortedM]=sort(m(:));
nv=numel(pctls);
vals=zeros(1,nv);
for i=1:nv
    ptr=round((np-1)*pctls(i)+1); % Get the index
    vals(i)=sortedM(max(1,min(ptr,numel(sortedM))));
end;
end


function h=NewAutomask(h,active)
if nargin<2
    active=true;
end;
if ~active
    %     h.mi=meInsertMask(true(round(size(h.origImage)/2)),h.mi,3,'AND');
    h.mi=InsertOurMask(true(size(h.origImage)/2),h.mi,h.ds0*2,3,'AND');
    h.miChanged=1;
    h=ShowImage(h);
    return
end;
if h.automaskOn
    %     thr1=h.sav.automaskPars.thresh*(h.automaskLims(2)-h.automaskLims(1))...
    %         +h.automaskLims(1);
    thr2=2*h.sav.automaskPars.width-2; % local mask threshold
    thr3=h.sav.automaskPars.dense;     % global mask threshold
    mskLocal=h.ifImageComp>thr2;
    mskGlobal= (h.e1Map>thr3) & (h.varMap<(1-h.sav.automaskPars.var));
    fc1=exp(-h.sav.automaskPars.edge*3); % local mask smoothness
    fc2=fc1/5;
    if fc1<.9
        mskLocal=GaussFilt((GaussFilt(mskLocal,fc1)>.1),fc2)>.9;
    end;
    fc1=exp(-h.sav.automaskPars.thresh*3)*.1; % global mask border
    mskGlobal=(GaussFiltDCT(mskGlobal,fc1)>.99);
    
    % The automask always goes into position 3.
    %     h.mi=meInsertMask(mskLocal & mskGlobal,h.mi,3,'AND');
    h.mi=InsertOurMask(mskLocal & mskGlobal,h.mi,h.ds0*2,3,'AND');
    h.miChanged=1;
    h=ShowImage(h);
end;
end


% --- Executes on slider movement.
function slider_Local_Callback(hObject, eventdata, h)
h.sav.automaskPars.width=get(h.slider_Local,'value');
h=NewAutomask(h);
guidata(hObject, h);
end

% --- Executes on slider movement.
function slider_Threshold_Callback(hObject, eventdata, h)
h.sav.automaskPars.thresh=get(h.slider_Threshold,'value');
h=NewAutomask(h);
guidata(hObject, h);
end

% --- Executes on slider movement.
function slider_GlobalVar_Callback(hObject, eventdata, h)
h.sav.automaskPars.var=get(h.slider_GlobalVar,'value');
h=NewAutomask(h);
guidata(hObject, h);
end

% --- Executes on slider movement.
function slider_Edges_Callback(hObject, eventdata, h)
h.sav.automaskPars.edge=get(h.slider_Edges,'value');
h=NewAutomask(h);
guidata(hObject, h);
end

% --- Executes on slider movement.
function slider_dense_Callback(hObject, eventdata, h)
h.sav.automaskPars.dense=get(h.slider_dense,'value');
h=NewAutomask(h);
guidata(hObject,h);
end




% --------AutomaskBeam---------

function h=UpdateAutomaskBeam(h,active)
% The beam mask always goes in index 2.
    h=ShowImage(h);
% 
% if nargin>1
%     h.automaskBeamOn=active;
% end;
% if h.automaskBeamOn
%     if ~isfield(h.mi,'mask')
%         h.mi.mask=struct('merge',[],'encoding',[],'data',[]);
%     end;
%     h.mi.mask(2).merge='AND';
%     h.mi.mask(2).encoding='RIM';
%     h.mi.mask(2).data=h.sav.beamPars/100;
%     h.miChanged=1;
%     %             return
%     %
%     %             %%%%%%%%%%
%     h=ShowImage(h);
% else
%     h.mi.mask(2).merge=[];
%     h=ShowImage(h);
% end;
end

function togglebutton_Beam_Callback(hObject, eventdata, h)
active=get(h.togglebutton_Beam,'value') && h.imageLoaded;
if active
    h.maskIndex=max(h.maskIndex,2);  % includes the beam at least.
end;
h=UpdateAutomaskBeam(h,active);
guidata(hObject, h);
end

function edit_beamX_Callback(hObject, eventdata, h)
q=str2double(get(hObject,'String'));
if ~isnan(q)
    h.sav.beamPars(1)=q;
    set(hObject,'string',num2str(q));
    h=UpdateAutomaskBeam(h);
    guidata(hObject, h);
end;
end

function edit_beamY_Callback(hObject, eventdata, h)
q=str2double(get(hObject,'String'));
if ~isnan(q)
    h.sav.beamPars(2)=q;
    set(hObject,'string',num2str(q));
    h=UpdateAutomaskBeam(h);
    guidata(hObject, h);
end;
end

function edit_beamR_Callback(hObject, eventdata, h)
q=str2double(get(hObject,'String'));
if ~isnan(q)
    h.sav.beamPars(3)=q;
    set(hObject,'string',num2str(q));
    h=UpdateAutomaskBeam(h);
    guidata(hObject, h);
end;
end



% hObject    handle to slider_dense (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%%%%%%%%%%%%%%%% FIND  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_Find.-----------------
function pushbutton_Find_Callback(hObject, eventdata, h)
h.doTrackMembranes=0;
DoFind(hObject,eventdata,h,false);
end

% --- Executes on button press in pushbutton_FindAndTrack.
function pushbutton_FindAndTrack_Callback(hObject, eventdata, h)
h.doTrackMembranes=1;
DoFind(hObject,eventdata,h,false);
end


% --- Executes on button press in pushbutton_FindMore.-----------------
% Searches the subtracted micrograph for more vesicles
function pushbutton_FindMore_Callback(hObject, eventdata, h)
h.doTrackMembranes=0;
DoFind(hObject,eventdata,h,true);
end

function DoFind(hObject, eventdata, h, findMore)
% Automatic finder. If findMore=true, it attempts to find vesicles in the
% image after vesicle models have been subtracted. The additional found
% vesicles are added to the mi.vesicle structure.
if ~h.imageLoaded
    return
end;
set(h.figure1,'pointer','watch');
% membrane model
% vLipid=h.sav.membranePars(1);
% thk=h.sav.membranePars(2);
% rise=h.sav.membranePars(3);
% pixA=h.mi.pixA;
%
% Load the generic model
% nm0=ceil(thk/(2*pixA))*2+3;  % array for vesicle model; 60A nominal
% if ~isfield(h.mi,'vesicleModel') || numel(h.mi.vesicleModel)<3
%     h.mi.vesicleModel=fuzzymask(nm0,1,thk/pixA/2,rise/pixA)...
%         *vLipid;  % units of V
% end;
% disp(' ');

% Instead, always load our standard membrane model.
pa=AddSlash(fileparts(which('Vesicle_finding_GUI'))); % stored with our code.
vm=load([pa 'VFDefaultVM.mat']);
h.mi.vesicleModel=meDownsampleVesicleModel(vm.vm1,h.mi.pixA);

rPars=h.sav.vesicleRadii;
mi1=h.mi;

if ~findMore
    % Zero out the previous particles if requested
    if h.sav.eraseOldPicks
        mi1.particle=struct;
        mi1.particle.picks=[];
        mi1.particle.autopickPars=[];
    end
    
    mi1.vesicle.x=[];
    mi1.vesicle.y=[];
    mi1.vesicle.r=[];
    mi1.vesicle.s=[];
    mi1.vesicle.ok=[];
    mi1.vesicle.err=[];
    if isfield(mi1,'mask')
        mi1.mask=mi1.mask(1:h.maskIndex); % collapse the mask array
    end;
    if ~isfield(mi1,'mergeMode')
        mi1.mergeMode=3;
    end;
end;
prevNFound=numel(mi1.vesicle.x);
pars.rPars=rPars;
pars.ds0=h.ds0;
pars.ds0Shift=-h.M2(1:2,3)';
if findMore
    mi1=rsFindVesicles3(h.origImage-h.rawVesImage, h.mi, pars, h.findInMask);
else
    mi1=rsFindVesicles3(h.origImage, mi1, pars, h.findInMask);
end;
%%
minAmp=h.sav.vesicleAmps(1);
maxAmp=h.sav.vesicleAmps(2);

mins=1;
nVesOld=0;
extraTries=2;
% Loop through finding groups of 50 vesicles.
while mins>minAmp
    [mi1, t]=rsFindVesicles3('next',50,h.sav.vesicleAmps);
    axes(h.axes1);
    ih = imshow(rot90(uint8(imscale(t.ms-t.umodel,256,1e-3))),...
        'InitialMagnification',100*size(t.ms,1)/size(h.origImage,1) );
    %     drawnow;
    mins=t.globalmax;
    nves=size(mi1.vesicle.s,1);
    set(hObject,'string',num2str(nves));
    %     Make the scatterplot of vesicle radii and amplitudes.
    
    %%     This scatterplot should be its own function...
    %     Make a rectangle to show the selected range of radius and amplitude values
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
    xlabel('Radius, A');
    ylabel('Amplitude');
    title(nves);
    %%
    
    %     Exit the loop when we can't find any more vesicles
    if nves<=nVesOld 
        extraTries=extraTries-1;
        disp('extra try');
        if extraTries<1 || nves > h.maxNVes
        
        break
        end;
    end;
    nVesOld=nves;
end;

% Store NCCs for manual finding
if ~findMore
    h.ccVals=0;
    h.ccValsScaled=0;
end;
h.ccVals=Crop(t.ccsmx,h.displaySize)+h.ccVals; % used for finding maxima
h.ccValsScaled=Crop(t.ccsmxScaled,h.displaySize)+h.ccValsScaled; % used for finding amplitude
h.ccRadii=(t.fitmin+(Crop(t.ccsmi,h.displaySize)-1)*t.rstep)*h.ds0;  % identify radius in orig pixels

mi1.vesicle.shiftX=[];
mi1.vesicle.shiftY=[];
mi1.vesicle.shiftOk=[];
mi1.vesicle.af=[];
mi1.vesicle.refined=0;
mi1.vesicle.extraPeaks=[];
mi1.vesicle.extraSD=0;
mi1.vesicle.extraS=[];

if findMore
    % Additional vesicles are marked 'false' in the 4th column.
    totalNFound=numel(mi1.vesicle.x);
    if totalNFound>prevNFound
        disp([num2str(totalNFound-prevNFound) ' additional vesicles found...']);
        mi1.vesicle.ok(prevNFound+1:totalNFound,4)=false;
    end;
end;

h.mi=mi1;
h.miChanged=1;
if h.makeModelVesicles
    scl.n=h.displaySize;
    scl.M=h.M2;
    h.goodVesImage=meMakeModelVesicles(h.mi,scl,find(goodVes));
    h.badVesImage=meMakeModelVesicles(h.mi,scl,find(badVes));
    h.rawVesImage=h.goodVesImage+h.badVesImage;
    % h.rawVesImage=meMakeModelVesicles(h.mi,size(h.origImage));
    h=UpdateDisplayFiltering(h);
end;
% h.displayMode=0;  % mark the vesicles
h.markedVesicleIndex=0;
ShowImage(h);
drawnow;

if h.doTrackMembranes && ~findMore
    
    % Tune up the vesicle fits
    disp('Tracking vesicle membranes');
    
    h.mi=TrackVesicleMembrane(h);
    
    % display everything again
    if h.makeModelVesicles
        disp('Computing vesicle models');
        goodVes=all(h.mi.vesicle.ok(:,1:2),2); % vesicles in range
        badVes=(h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2)); % found, but not in range
        h.goodVesImage=meMakeModelVesicles(h.mi,scl,find(goodVes));
        h.badVesImage=meMakeModelVesicles(h.mi,scl,find(badVes));
        h.rawVesImage=h.goodVesImage+h.badVesImage;
        % h.rawVesImage=meMakeModelVesicles(h.mi,size(h.origImage));
        h=UpdateDisplayFiltering(h);
    end;
    ShowImage(h);
    drawnow;
    
end;

guidata(hObject,h);
set(hObject,'string','Find');
set(h.figure1,'pointer','arrow');
disp('Done finding.');
end


% function image_ButtonDownFcn(hObjec,~,handles);
% % hObject    handle to figure1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% disp('image')
%     p=get(handles.axes1,'CurrentPoint')
% end
%

%


%%%%%%%%%%%%% MOUSE BUTTONS %%%%%%%%%%%%%%%%%%%%%%%%%

function axes1_ButtonDownFcn(hObject, eventdata, h)
axesHandle=get(hObject,'parent');
h=guidata(axesHandle);
q=get(h.pushbutton_OutlineMask,'value');

if h.manualMaskActive      % draw on the screen
    if h.manualMaskDiameter>0
        set(h.axes1,'nextplot','add');
        pt=get(h.axes1,'currentpoint');
        h.manualMaskCoords(end+1,:)=pt(1,1:2);
        plot(h.axes1,pt(1,1),pt(1,2),'o','markersize',h.manualMaskDiameter,'markerfacecolor','b');
        set(h.axes1,'nextplot','add');
        h.axes1ButtonDown=1;
    end;
elseif ~q
    [h, doUpdateDisplay]=vfMouseClick(h);
    if doUpdateDisplay
        h=UpdateDisplayFiltering(h,false);
    end;
    h=ShowImage(h);
end
guidata(axesHandle,h);
end


function WindowButtonMotionFcn(hObject, eventdata)
h=guidata(hObject);
if h.axes1ButtonDown && h.manualMaskActive && h.manualMaskDiameter>0
    pt=get(h.axes1,'currentpoint');
    h.manualMaskCoords(end+1,:)=pt(1,1:2);
    %         if numel(h.theImage)>0
    %             n=size(h.theImage);
    %             msk=(1-uint8(disc(n,50,pt(1,:))))';
    %             h.theImage(:,:,1)=h.theImage(:,:,1).*msk;
    %                 axes(h.axes1);
    %             imshow(h.theImage,'InitialMagnification',100);
    %         end;
    set(h.axes1,'nextplot','add');
    plot(h.axes1,pt(1,1),pt(1,2),'o','markersize',h.manualMaskDiameter,'markerfacecolor','b');
    set(h.axes1,'nextplot','add');
    
    %         pt=get(h.axes3,'currentpoint');
    %         disp(pt);
    guidata(hObject,h);
end
end

function WindowButtonUpFcn(hObject, eventdata)
h=guidata(hObject);

if h.manualMaskActive && h.manualMaskDiameter>0
    h=InsertPaintedMask(h);
end;

h.axes1ButtonDown=false;
guidata(hObject,h);
% pt=get(hObject,'currentpoint');
% disp(pt);
end


%---------------keystroke------------------------
% For some reason, user has to click on the figure first for this to work.
function KeyPressFcn(hObject, eventdata, h)
figure(gcbf);
char=get(gcbf,'CurrentCharacter');
switch char
    case {'d' ' '}
        h.displayMode=h.displayMode+1;
        if h.displayMode>h.maxDisplayMode
            h.displayMode=0;
        end;
        h=ShowImage(h);
end;
guidata(hObject,h);

end



% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, h)
pos=get(hObject,'currentpoint');
disp(pos);

% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
% Vesicle_finding_GUI('axes3_ButtonDownFcn',hObject,eventdata,guidata(hObject))

% --- Executes on button press in pushbutton_Save.
function pushbutton_Save_Callback(hObject, eventdata, h)

end

% --- Executes on button press in pushbutton_Revert.
function pushbutton_Revert_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton_Revert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
end


function edit_Edge_Callback(hObject, eventdata, h)
val=str2double(get(hObject,'String')); % returns contents of edit_Edge as a double
if val>0 && val<40
    h.sav.membranePars(3)=val;
    guidata(hObject,h);
end;
end

function edit_Thickness_Callback(hObject, eventdata, h)
val=str2double(get(hObject,'String')); % returns contents of edit_Edge as a double
if val>=10 && val<=100
    h.sav.membranePars(2)=val;
    guidata(hObject,h);
end;
end

% --- Executes on button press in pushbutton_LoadModel.
function pushbutton_LoadModel_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton_LoadModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function edit_MinAmp_Callback(hObject, eventdata, h)
val=str2double(get(hObject,'String')); % returns contents of edit_MinRadius as a double
if val > 0 && val < 1000
    h.sav.vesicleAmps(1)=val*1e-3;
    guidata(hObject,h);
end;
end

function edit_MaxAmp_Callback(hObject, eventdata, h)
val=str2double(get(hObject,'String')); % returns contents of edit_MinRadius as a double
if val > 0 && val < 1000
    h.sav.vesicleAmps(2)=val*1e-3;
    guidata(hObject,h);
end;
end

% We now use iceAmp as the vesicle boost radius. We don't scale the value.
function edit_IceAmp_Callback(hObject, eventdata, h)
val=str2double(get(hObject,'String')); % returns contents of edit_IceAmp as a double
if val >= 0 && val <= 1
    h.sav.vesicleAmps(3)=val;
    guidata(hObject,h);
end;
end


function edit_MinRadius_Callback(hObject, eventdata, h)
val=str2double(get(hObject,'String')); % returns contents of edit_MinRadius as a double
if val > 20 && val < 1000
    h.sav.vesicleRadii(1)=val;
    guidata(hObject,h);
end;
end

function edit_MaxRadius_Callback(hObject, eventdata, h)
val=str2double(get(hObject,'String')); % returns contents of edit_MinRadius as a double
if val > 100 && val < 2000
    h.sav.vesicleRadii(2)=val;
    guidata(hObject,h);
end;
end

function edit_RadiusStep_Callback(hObject, eventdata, h)
val=str2double(get(hObject,'String')); % returns contents of edit_MinRadius as a double
if val > 1 && val < 100
    h.sav.vesicleRadii(3)=val;
    guidata(hObject,h);
end;
end

% --- Executes on button press in pushbutton_ShowVesicles.
function pushbutton_ShowVesicles_Callback(hObject, eventdata, h)
h.displayMode=h.displayMode+1;
if h.displayMode>h.maxDisplayMode
    h.displayMode=0;
end;
h=ShowImage(h);
guidata(hObject,h);
end


% --- Executes on button press in FindInMaskButton.
function FindInMaskButton_Callback(hObject, eventdata, h)
% hObject    handle to FindInMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h.findInMask=1-h.findInMask;
set(hObject,'Value',h.findInMask);
guidata(hObject,h);
end
% Hint: get(hObject,'Value') returns toggle state of FindInMaskButton


% --- Executes on button press in EraseOldPicksButton.
function EraseOldPicksButton_Callback(hObject, eventdata, h)
% hObject    handle to EraseOldPicksButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h.sav.eraseOldPicks=1-h.sav.eraseOldPicks;
set(hObject,'Value',h.sav.eraseOldPicks);
guidata(hObject,h);
end


% --- Executes on key press with focus on pushbutton_AdjustContrast and none of its controls.
function pushbutton_AdjustContrast_KeyPressFcn(hObject, eventdata, h)
% hObject    handle to pushbutton_AdjustContrast (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% h    structure with h and user data (see GUIDATA)
end

% --- Executes on button press in pushbutton_DeleteAll.
function pushbutton_DeleteAll_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton_DeleteAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
end

% --- Executes on button press in pushbutton_Delete.
function pushbutton_Delete_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton_Delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
end

% --- Executes on button press in pushbutton_select.
function pushbutton_Select_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
end

% --- Executes on key press with focus on pushbutton_select and none of its controls.
function pushbutton_Select_KeyPressFcn(hObject, eventdata, h)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% h    structure with h and user data (see GUIDATA)
end

% --- Executes on key press with focus on pushbutton_Find and none of its controls.
function pushbutton_Find_KeyPressFcn(hObject, eventdata, h)
% hObject    handle to pushbutton_Find (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% h    structure with h and user data (see GUIDATA)
end




function FileNameText_Callback(hObject, eventdata, h)
% hObject    handle to FileNameText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileNameText as text
%        str2double(get(hObject,'String')) returns contents of FileNameText as a double

end
% --- Executes during object creation, after setting all properties.
function FileNameText_CreateFcn(hObject, eventdata, h)
% hObject    handle to FileNameText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



% --- Executes during object creation, after setting all properties.
function slider_Local_CreateFcn(hObject, eventdata, h)
% hObject    handle to slider_Local (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end

% --- Executes during object creation, after setting all properties.
function slider_Threshold_CreateFcn(hObject, eventdata, h)
% hObject    handle to slider_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes during object creation, after setting all properties.
function slider_Edges_CreateFcn(hObject, eventdata, h)
% hObject    handle to slider_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes when Main_im is resized.
function Main_im_ResizeFcn(hObject, eventdata, h)
% hObject    handle to Main_im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
end



% --- Executes during object creation, after setting all properties.
function edit_Highpass_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_Highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function edit_Lowpass_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_Lowpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Hand-written callback
% --- Used to return 'CData' for the Stop icon on the Record\Stop toggle button
function invarrow = inarrowbutton

invarrow = iconize(imread('inverse_arrow.jpg'));
invarrow(invarrow==255) = .8*255;
end

% --- Hand-written callback
% --- Used to return 'CData' for the Play icon on the Play button
function arrow = arrowbutton

arrow = iconize(imread('arrow.jpg'));
arrow(arrow==255) = .8*255;
end

% --- Hand-written callback
% --- Used to create icon data from an image, a
function out = iconize(a)
% Find the size of the acquired image and determine how much data will need
% to be lost in order to form a 18x18 icon
[r,c,d] = size(a);
r_skip = ceil(r/9);
c_skip = ceil(c/9);
% Create the 18x18 icon (RGB data)
out =  a(1:r_skip:end,1:c_skip:end,:);
end

% --- Executes during object creation, after setting all properties.
function edit_Edge_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_Edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function edit_MinRadius_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_MinRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes during object creation, after setting all properties.
function edit_MaxRadius_CreateFcn(hObject, ~, h)
% hObject    handle to edit_MaxRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes during object creation, after setting all properties.
function edit_RadiusStep_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_RadiusStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, h)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
end

% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, h)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, h)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
end

% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, h)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end



% --- Executes during object creation, after setting all properties.
function edit_MinAmp_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_MinAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function edit_MaxAmp_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_MaxAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end




% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
disp('pushbutton22');
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes during object creation, after setting all properties.
function slider_dense_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_dense (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes during object creation, after setting all properties.
function edit_beamX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beamX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function edit_beamY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beamY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function edit_beamR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beamR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
%
% end


% --- Executes during object creation, after setting all properties.
function editWhite_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWhite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;
end

% --- Executes during object creation, after setting all properties.
function editBlack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBlack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;
end


% --- Executes during object creation, after setting all properties.
function edit_IceAmp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IceAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function slider_GlobalVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_GlobalVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end
