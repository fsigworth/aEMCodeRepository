function varargout = Vesicle_finding_GUI(varargin)
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

% Last Modified by GUIDE v2.5 29-Jul-2013 19:16:54

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
screenSize = get(0, 'ScreenSize');

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
amPars.thresh=.5;
amPars.width=.5;
amPars.edge=.5;
sav.automaskPars=amPars;  % saved variables, stored in "VFContext.mat"
sav.basePath='';
sav.fullInfoPath='';
sav.baseName='';
sav.vesicleAmps=[5e-3 1e-3];
sav.vesicleRadii=[100 300 10];
sav.filterFreqs=[.01 .1];   % HP and LP, in A^-1
sav.contrastPars=zeros(2,1);
sav.membranePars=[1.6 60 6];
% State variables
h.sav=sav;
h.imageLoaded=false;
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
h.oldFilterFreqs=[0 0];
h.maskIndex=3;       % default uses existing masks.
h.vesModels=single(0);
h.ds=1;
h.miChanged=false;
h.mi=struct;
h.miOriginal=struct;
h.automaskOn=false;
h.varianceMap=single(0);
h.axes1ButtonDown=false;
h.markedVesicleIndex=0;
h.findInMask=1;
h.eraseOldPicks=1;
h.borderFraction=256;  % relative size of merged-image border is 1/256 ------------


pos=get(h.axes1,'position');
h.displaySize=pos(3:4);
h.displaySize=768;

% h for the image display
h.hSP=0;
h.hNav=0;
h.ih=0;
% h.sav.automaskPars

%%% read a context file to get the parameters
% Find our local directory
pa=fileparts(which('Vesicle_finding_GUI'));
% Retrieve parameters from VesiclePara.mat in the local directory
if exist([pa '/VFContext.mat'],'file')
    disp('Loading VFContext.mat');
    sav=load([pa '/VFContext.mat']);
    h.sav=sav.sav;
end;

% set the slider parameters
set(h.slider_Threshold,'Min', 0);
set(h.slider_Threshold,'Max', 1.5);
set(h.slider_Threshold,'SliderStep',[1/(200-1) .2]);
set(h.slider_Threshold,'Value',h.sav.automaskPars.thresh);

set(h.slider_Smoothness,'Min', 0);
set(h.slider_Smoothness,'Max', 1);
set(h.slider_Smoothness,'SliderStep',[1/(200-1) .2]);
set(h.slider_Smoothness,'Value',h.sav.automaskPars.width);

set(h.slider_Edges,'Min', 0);
set(h.slider_Edges,'Max', 1);
set(h.slider_Edges,'SliderStep',[1/(200-1) .2]);
set(h.slider_Edges,'Value',h.sav.automaskPars.edge);

% set the membrane parameters
% set(h.edit_Edge,'String',sprintf('%1.1f',h.sav.membranePars(3)));
% set(h.edit_Thickness,'String',sprintf('%d',h.sav.membranePars(2)));

set(h.edit_MinRadius,'String',sprintf('%d',h.sav.vesicleRadii(1)));
set(h.edit_MaxRadius,'String',sprintf('%d',h.sav.vesicleRadii(2)));
set(h.edit_MinAmp,'String',num2str(1000*h.sav.vesicleAmps(1)));
set(h.edit_MaxAmp,'String',num2str(1000*h.sav.vesicleAmps(2)));
set(h.edit_RadiusStep,'String',sprintf('%d',h.sav.vesicleRadii(3)));
set(h.edit_Highpass,'string',num2str(h.sav.filterFreqs(1)));
set(h.edit_Lowpass,'string',num2str(h.sav.filterFreqs(2)));

set(h.FindInMaskButton,'value',h.findInMask);
set(h.EraseOldPicksButton,'value',h.findInMask);

h.output = hObject;
h.altBasePath='';  % mi.basePath value if the original one doesn't work.
h.imageFileTypes={'m.mrc' 'mf.mrc' 'm.jpg' 'mf.jpg'};  % allowed types of
% images to load based on mi file.

% h
% Update h structure
guidata(hObject, h);
% UIWAIT makes Vesicle_finding_GUI wait for user response (see UIRESUME)
% uiwait(h.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = Vesicle_finding_GUI_OutputFcn(~, ~, h)
varargout{1} = h.output;
end


function figure1_CloseRequestFcn(hObject, eventdata, h)
pa=fileparts(which('Vesicle_finding_GUI'));
% Save parameters from VesiclePara.mat in the local directory
sav=h.sav;
save([pa '/VFContext.mat'],'sav');
CloseFile(h);
rsFindVesicles3('end');  % purge the variables
delete(hObject);
end


function CloseFile(h) % store the results of operations
if ~h.imageLoaded
    return
end;
% discard unused masks
if isfield(h.mi,'mask')&&(numel(h.mi.mask)>h.maskIndex)
    h.mi.mask=h.mi.mask(1:h.maskIndex);
end;
% mi=h.mi;
if h.eraseOldPicks
    h.mi.particle.picks=[];
    mi=h.mi;
else
    mi=rsMergeVesicleList2(h.mi,h.miOriginal);
end;
fileName=[h.sav.baseName 'mi.mat'];
outName=[h.sav.fullInfoPath fileName];
disp(['Saving: ' outName]);
save(outName,'mi');  % stores mi
end


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% --- Executes on button press in pushbutton_LoadFile.
function pushbutton_LoadFile_Callback(hObject, ~, h)
% We'll read an mi file
if ~isempty(h.sav.fullInfoPath)&& exist(h.sav.fullInfoPath,'dir')
    disp(['full info path is ' h.sav.fullInfoPath]);
    tempPath=pwd;
    cd(h.sav.fullInfoPath);
    [fileName,infoPath] = uigetfile({'*.mat'},'Load Data File');
    cd(tempPath);
else
    [fileName,infoPath] = uigetfile({'*.mat' },'Load Data File');
end
if (fileName==0)
    return
end;
[~, nm, ex] = fileparts(fileName);
len=numel(nm);
if  len<2 || ~strcmpi(ex,'.mat') || ~strcmp(nm(len-1:len),'mi')   % judge whether choose mat file
    msgbox('Expected mi.mat file.','Error','error');
    return;
end

CloseFile(h);  % Save the previous file.

h.sav.fullInfoPath=AddSlash(infoPath);
mi=load([h.sav.fullInfoPath fileName]);  % loads mi
mi=mi.mi;
set(h.FileNameText,'String',fileName);
h.sav.baseName=nm(1:numel(nm)-2);  % delete the 'mi'


[h, ok]=GetImageFile(mi,h);
if ~ok
    h.imageLoaded=false;
    return
end;

h=InitDisplay(h);
guidata(hObject, h);
end



function [h ok]=GetImageFile(mi,h)
% Look for a merged image file.


% update the vesicle info to the canonical nv x 4 logical matrix
nv=numel(mi.vesicle.x);
if (~isfield(mi.vesicle,'ok')...
        && size(mi.vesicle.ok,1)==nv)
    mi.vesicle.ok=true(nv,1);
end;
if size(mi.vesicle.ok,2)<4 && nv>0
    mi.vesicle.ok(nv,4)=false;  % extend the array
end;


basePath=mi.basePath;
if ~exist(basePath,'dir')  % Original directory doesn't work, try current dir
    mi.basePath=ParsePath(h.sav.fullInfoPath);
    if ~exist(mi.basePath,'dir')  % Can't find base directory
        msgbox(['Can''t find the path ' mi.basePath '.   Please pick the base directory'],'ok');
        mi.basePath=AddSlash(uigetdir(ParsePath(h.currentInfoPath)));
        if isnumeric(basePath)  % cancel was clicked
            ok=0;
            return
        end;
    end;
end;
%     Try for a merged file
imageBasename=[mi.basePath mi.procPath mi.baseFilename];
ok=0;
for i=1:numel(h.imageFileTypes)
    fullImageName=[imageBasename h.imageFileTypes{i}];
    if exist(fullImageName,'file')
        ok=1;
        break;
    end;
end;
if ok  % valid filename.  Load it
    h.oldFilterFreqs=[0 0];
    h.sav.basePath=mi.basePath;
    h.baseName=imageBasename;
    %     h.origImage=ReadEMFile(fullImageName)/mi.doses(1);  % normalize by dose.
    h.origImage=ReadEMFile(fullImageName);
    h.imageLoaded=true;
    h.rawImage=single(0);
    h.automaskOn=false;
    set(h.togglebutton_Automask,'value',0);
    %     Clear all the processed images.
    h.filtImage=single(0);
    h.rawVesImage=single(0);
    h.filtVesImage=single(0);
    cla(h.axes3);
% initialize the mask
    t=mi.mergeMatrix;
    msk=meMakeMergedImageMask(mi.imageSize/4,t,mi.imageSize/(4*h.borderFraction));
    mi=meInsertMask(msk,mi,1);
%     Point to the end of the stack
h.maskIndex=numel(mi.mask);

%   Initialize variables
h.miOriginal=mi;
h.mi=mi;
h.sav.basePath=mi.basePath;
h.ccVals=0;
h.ccValsScaled=0;
h.ccRadii=0;
h.markedVesicleIndex=0;

else
    msgbox(['Can''t find the image ' imageBasename],'ok');
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
disp('Making model vesicles')

[nv, ne]=size(h.mi.vesicle.ok);
if ne<4
    h.mi.vesicle.ok(1,4)=false;  % extend it.
end;

goodVes=all(h.mi.vesicle.ok(:,1:2),2); % vesicles in range
badVes=(h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2)); % found, but not in range

h.goodVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(goodVes));
h.badVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(badVes));
h.ctf=meGetEffectiveCTF(h.mi,size(h.rawImage),h.ds0);
h=UpdateDisplayFiltering(h);
% h.displayMode=0;
h=ShowImage(h);
% set(h.figure1,'pointer','arrow');

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


function h=UpdateDisplayMasking(h)
h.rootMask=meGetMask(h.mi,h.displaySize);
end


% ----ShowImage----
function h=ShowImage(h)
if h.imageLoaded
    % if numel(h.filtImage)>0  % something has been loaded
    msk=rot90(meGetMask(h.mi,size(h.filtImage),1:h.maskIndex));
    imData=h.filtImage;
    n=size(h.filtImage);
    showMask=1;
    switch h.displayMode
        case 0
            imData=h.filtImage;
            showGhosts=0;
        case 1
            imData=h.filtImage-h.filtVesImage;
            showGhosts=0;
        case 2
            imData=h.filtImage-h.filtVesImage;
            showGhosts=1;
        case 3
            if numel(h.ccValsScaled)>1
                imData=h.ccValsScaled(1:n(1),1:n(2));
                showGhosts=1;
            else
                return
            end;
    end;
    theImage =  repmat(rot90(imscale(imData,256,1e-3)),[1 1 3]);
    
    if showGhosts
        ghostColor=[.7 .7 1];
        ghostColorBad=[1 .5 .35];
        
        gves=imscale(max(-h.goodVesImage,0),1);
        bves=imscale(max(-h.badVesImage,0),1);
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
        maskColor=[.7 .6 .7];
        
        color=(1-maskColor);
        for i=1:3
            theImage(:,:,i)=theImage(:,:,i).*(1-(1-msk)*color(i));
        end;
    end;
    axes(h.axes1);
    ih = imshow(uint8(theImage),'InitialMagnification',100);
    if showGhosts % draw the center points of the vesicles
        %         goodCenterColor=[.5 .5 1];
        %         badCenterColor=[1 .5 1];
        if isfield(h.mi.vesicle,'ok') && size(h.mi.vesicle.ok,2)>3
            badVes=h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2);
            goodVes=all(h.mi.vesicle.ok(:,1:2),2);
        else
            badVes=false(numel(h.mi.vesicle.x),1);
            goodVes=true(numel(h.mi.vesicle.x),1);
        end;
        hold on;
        ny=size(h.rawImage,2);
        hp1=plot(h.mi.vesicle.x(goodVes)/h.ds0+1,...
            ny-(h.mi.vesicle.y(goodVes)/h.ds0),'b.','markersize',10);
        hp2=plot(h.mi.vesicle.x(badVes)/h.ds0+1,...
            ny-(h.mi.vesicle.y(badVes)/h.ds0),'r.','markersize',10);
        hold off;
        set(hp1,'HitTest','off');
        set(hp2,'HitTest','off');
    end;
    
    % Draw the scatterplot of existing vesicles
    
    if isfield(h.mi.vesicle,'ok')
        [nv, ne]=size(h.mi.vesicle.ok);
        if nv>0 && ne>3
            goodVes=all(h.mi.vesicle.ok(:,1:2),2); % vesicles in range
            badVes=(h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2)); % found, but not in range
            plot(h.axes3,...
                h.mi.vesicle.r(goodVes)*h.mi.pixA,h.mi.vesicle.s(goodVes),'b.',...
                h.mi.vesicle.r(badVes)*h.mi.pixA,h.mi.vesicle.s(badVes),'m.',...
                'markersize',10);
            if h.markedVesicleIndex>0
                set(h.axes3,'nextplot','add');
                plot(h.axes3,h.mi.vesicle.r(h.markedVesicleIndex)*h.mi.pixA,...
                    h.mi.vesicle.s(h.markedVesicleIndex),'ks','markersize',12);
                set(h.axes3,'nextplot','replace');
            end;
            title(h.axes3,[num2str(sum(goodVes)) ' vesicles in range']);
        else
            title(h.axes3,' '); % clear the title
        end;
    end;
    
    fhndl=@(hObject,eventdata) Vesicle_finding_GUI('axes1_ButtonDownFcn',hObject,eventdata,guidata(hObject));
    set(ih,'buttondownfcn',fhndl);
    set(ih,'HitTest','on');
    h.ih=ih;
end
end

% Vesicle_finding_GUI('axes1_ButtonDownFcn',hObject,eventdata,guidata(hObject))


% --- Executes on button press in pushbutton_formerName.
function pushbutton_formerName_Callback(hObject, ~, h)
[h ok]=LoadAnotherMiFile(h,-1);
if ok
    guidata(hObject, h);
end;
end


% --- Executes on button press in pushbutton_nextName.
function pushbutton_nextName_Callback(hObject, ~, h)
[h ok]=LoadAnotherMiFile(h,1);
if ok
    guidata(hObject, h);
end;
end




function [h ok]=LoadAnotherMiFile(h,offset)
ok=0;
if ~h.imageLoaded
    return
end;
d=dir(h.sav.fullInfoPath);
dirIndex = [d.isdir];  %# Find the index for directories
fileList = {d(~dirIndex).name}';  %'# Get a list of the files

findcurrentfile = strfind(fileList,[h.sav.baseName 'mi.mat']);
currentfile_id = find(~cellfun(@isempty,findcurrentfile));   % current working file id
id=currentfile_id+offset;
if id<1 || id>numel(fileList)
    beep;
    return;
end
fileName=fileList{id};
CloseFile(h);  % save the previous one.
mi=load([h.sav.fullInfoPath fileName]) ;
set(h.FileNameText,'String',fileName);

mi=mi.mi;
[pa nm ex]=fileparts(fileName);
h.sav.baseName=nm(1:numel(nm)-2);  % delete the 'mi'

[h ok]=GetImageFile(mi,h);
if ~ok
    return
end;

% Update the mask

% if ~isfield(mi,'mask') || ~isfield(mi.mask,'merge') || numel(mi.mask.merge)<1
%    Old merged data: compute the base mask and insert it.
% disp('Computing merge mask...');
t=mi.mergeMatrix;
msk=meMakeMergedImageMask(mi.imageSize/4,t,mi.imageSize/(4*h.borderFraction));
mi=meInsertMask(msk,mi,1);
% end;
%     Point to the end of the stack
h.maskIndex=numel(mi.mask);

h.mi=mi;  % commit to processing this image.
h.miOriginal=mi;
h.sav.basePath=mi.basePath;

h=InitDisplay(h);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_ManualMask.
function pushbutton_ManualMask_Callback(hObject, eventdata, h)
% Use the stored property 'value' in the axes1 button down function to
% actually initiate the manual mask input.
if ~h.imageLoaded
    set(h.pushbutton_ManualMask,'value',0);
elseif ~isfield(h.mi,'mask')
    h.mi.mask=struct('merge',[],'encoding',[],'data',[]);
end;
q=get(h.pushbutton_ManualMask,'value');
if h.imageLoaded && q % We're doing manual masking
    try
        hFH = imfreehand();
        binaryImage =~rot90(hFH.createMask(),3);
        h.maskIndex=max(3,h.maskIndex)+1;
        h.mi=meInsertMask(binaryImage,h.mi,h.maskIndex,'AND');
        h=ShowImage(h);
    catch   % Error occurred, exit manual mask mode.
        set(h.pushbutton_ManualMask,'value',0);
        return
    end;
    
    % Create a binary image ("mask") from the ROI object.
    set(h.pushbutton_ManualMask,'value',0);
    guidata(hObject, h);
end;
end


% --- Executes on button press in pushbutton_Invert.
function pushbutton_Invert_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton_Invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
msk=meGetMask(h.mi,0,h.maskIndex);  % get mask at native size.
h.mi=meInsertMask(~msk,h.mi,h.maskIndex,...
    h.mi.mask(h.maskIndex).merge);
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------Automask------------

% --- Executes on button press in togglebutton_Automask.
function togglebutton_Automask_Callback(hObject, eventdata, h)
h.automaskOn=get(h.togglebutton_Automask,'value') && h.imageLoaded;
if h.automaskOn
    if ~isfield(h.mi,'mask')
        h.mi.mask=struct('merge',[],'encoding',[],'data',[]);
    end;
    oPixA=h.mi.pixA*h.mi.imageSize(1)/size(h.origImage,1);
    varMode=get(h.checkbox_Var,'value');
    [msk h.varianceMap]=AutoHoleMask(h.origImage,oPixA,...
        h.sav.automaskPars,varMode);
    h=NewAutomask(h);
end;
guidata(hObject, h);
end

function h=NewAutomask(h)
if h.automaskOn
    h.maskIndex=max(h.maskIndex,3);  % we'll write mask 3
    oPixA=h.mi.pixA*h.mi.imageSize(1)/size(h.origImage,1);
    varMode=get(h.checkbox_Var,'value');
    msk=AutoHoleMask(h.origImage,oPixA,h.sav.automaskPars,...
        varMode,h.varianceMap);
    % The automask always goes into position 3.
    h.mi=meInsertMask(msk,h.mi,3,'AND');
    h=ShowImage(h);
end;
end

% --- Executes on button press in checkbox_Var.
function checkbox_Var_Callback(hObject, eventdata, h)

% hObject    handle to checkbox_Var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Var
end

% --- Executes on slider movement.
function slider_Smoothness_Callback(hObject, eventdata, h)
h.sav.automaskPars.width=get(h.slider_Smoothness,'value');
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
function slider_Edges_Callback(hObject, eventdata, h)
h.sav.automaskPars.edge=get(h.slider_Edges,'value');
h=NewAutomask(h);
guidata(hObject, h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_Find.-----------------
function pushbutton_Find_Callback(hObject, eventdata, h)
if ~h.imageLoaded
    return
end;

set(h.figure1,'pointer','watch');
% membrane model
vLipid=1.6;
thk=h.sav.membranePars(2);
rise=h.sav.membranePars(3);
pixA=h.mi.pixA;
% Create the model, which is sampled in units of the original pixel size.
nm0=ceil(30/pixA)*2+1;  % array for vesicle model; 60A nominal
h.mi.vesicleModel=fuzzymask(nm0,1,thk/pixA/2,rise/pixA)...
    *vLipid;  % units of V

rPars=h.sav.vesicleRadii;
% mi1=rsFindVesicles3(h.origImage, h.mi, rPars);
mi1=rsFindVesicles3(h.rawImage, h.mi, rPars, h.findInMask);
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
    nves=numel(mi1.vesicle.s);
    set(hObject,'string',num2str(nves));
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
    plot(h.axes3,mi1.vesicle.r(goodVes)*mi1.pixA,mi1.vesicle.s(goodVes),'b.',...
        mi1.vesicle.r(badVes)*mi1.pixA,mi1.vesicle.s(badVes),'m.',...
        xs,ys,'k-','markersize',10);
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
mi1.vesicle.refined=0;
h.mi=mi1;
h.goodVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(goodVes));
h.badVesImage=meMakeModelVesicles(h.mi,size(h.rawImage),find(badVes));
h.rawVesImage=h.goodVesImage+h.badVesImage;
% h.rawVesImage=meMakeModelVesicles(h.mi,size(h.rawImage));
h=UpdateDisplayFiltering(h);
h.displayMode=2;  % subtract and mark the vesicles
h.markedVesicleIndex=0;
ShowImage(h);
guidata(hObject,h);
set(hObject,'string','Find');
set(h.figure1,'pointer','arrow');
disp('Done.');
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

%---------------click on display-----------------

function axes1_ButtonDownFcn(hObject, eventdata, h)

axesHandle=get(hObject,'parent');
h=guidata(axesHandle);
[h, doUpdateDisplay]=vfMouseClick(h);
if doUpdateDisplay
    h=UpdateDisplayFiltering(h,false);
end;
h=ShowImage(h);
guidata(axesHandle,h);
end

%---------------keystroke------------------------
% For some reason, user has to click on the figure first for this to work.
function KeyPressFcn(hObject, eventdata, h)
figure(gcbf);
char=get(gcbf,'CurrentCharacter');
switch char
    case {'d' ' '}
        h.displayMode=h.displayMode+1;
        if h.displayMode>2
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
if h.displayMode>2
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
h.eraseOldPicks=1-h.eraseOldPicks;
set(hObject,'Value',h.eraseOldPicks);
guidata(hObject,h);
end

% Hint: get(hObject,'Value') returns toggle state of EraseOldPicksButton







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
function slider_Smoothness_CreateFcn(hObject, eventdata, h)
% hObject    handle to slider_Smoothness (see GCBO)
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
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

