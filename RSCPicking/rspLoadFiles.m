function [dis, mi, rscc, rs, imgs, masks, rawMask,success]=rspLoadFiles(dis,si)
% Given the info file name and path (in the dis structure), load the mi, merged image and rscc
% files.  If a vesicle file is also present load it, otherwise compute the
% vesicles.
dsc=4; % downsample factor for fake rscc sie
forceModelVesicles=0; % force the computing of vesicle models.

dis.basePath=AddSlash(pwd);
disp(' ');
if dis.miIndex
    disp(['file index: ' num2str(dis.miIndex)]);
end;
disp([dis.infoPath dis.infoName]);

if dis.miMode
    mi=ReadMiFile([AddSlash(dis.infoPath) dis.infoName],1);  % get the mi structure
else
    mi=si.mi{dis.miIndex};
end;
mi.basePath=dis.basePath;
set(gcf,'name',[dis.infoName ' - loading...']);
dis.miValid=1;

% handle bug where weights are zero.
if sum(mi.weights)==0
    mi.weights=ones(numel(mi.ctf),1);
end;

% Default returned values
rscc=struct;
rs=struct;
imgs=[];
masks=[];
rawMask=[];
success=false;

%sufExts={'s.mrc' 'z.tif' '.mrc'}; % suffix and extension options for small, compressed or full files.
% drawnow;
[mc, mergedName,ok]=meReadMergedImage(mi,0,'s');
%[mergedName,ok]=CheckForAltImage([mi.procPath mi.baseFilename 'm.mrc'],sufExts);
% [mc,mergedName]=meReadMergedImage(mi);
% mc=ReadEMFile(mergedName);
disp(mergedName);
if ~ok  % unsuccessful read
    mc=single(zeros(960));
    disp('No merged image found.');
    return
end;
if isfield(mi,'procPath_sm')
    procPath_sm=mi.procPath_sm;
else
    procPath_sm=mi.procPath;
end;
rsccName=[procPath_sm mi.baseFilename 'rscc.mat'];
disp(rsccName);
dis.rsMode=exist(rsccName,'file');
if dis.rsMode
    rscc=load(rsccName);
    % ---------check and set the window size-------
    nRscc=size(rscc.mxCC,1);
    if nRscc ~= dis.ndis(1)
        dis=ChangeTheWindowSize(dis,rscc);
    end;
    if isfield(rscc,'ppVals')
        mi.ppVals=rscc.ppVals;   % copy the parameters from the picking preprocessor
    end;
    if isfield(rscc,'spMode')
        dis.spMode=rscc.spMode;
    end
    if ~isfield(rscc,'mxCCU')
        rscc.mxCCU=rscc.mxCC*.001;
        disp('mxCCU not present, using mxCC.');
    end;
    
else % set up for generic particle picking
    disp('...not found');
    %     return
    
    
    
    if ~all(mi.imageSize>0)
        mi.imageSize=dis.ndis(1);
    end;
    rscc.spMode=1;
    rscc.mxVars=zeros(dis.ndis,'single');
    rscc.mxCC=rscc.mxVars;
    rscc.mxCCU=rscc.mxVars;
    rscc.mxVesInds=rscc.mxVars;
    rscc.mxTemplInds=rscc.mxVars;
    rscc.mxRsos=rscc.mxVars;
    rscc.mVes=rscc.mxVars;
    rscc.eigenImgs=zeros(32,32,20,'single');
end;

% initialize the CC for autopicking
if dis.useRawAmplitudes
    rscc.mxCC2=rscc.mxCCU*dis.ccuScale;
else
    rscc.mxCC2=rscc.mxCC;
end;

if ~isfield(rscc,'membraneOffsetA')
    rscc.membraneOffsetA=0;
end;
if ~isfield(mi,'particle')
    mi.particle=struct;
end;
if ~isfield(mi.particle,'picks')
    mi.particle.picks=[];
end;
if isfield(rscc,'pwfRef') % The prewhitening filter for particles
    ne=size(rscc.pwfRef,1);
    ct=ceil((ne+1)/2);
    rscc.pwfRef(ct,ct)=inf;
    dis.invPWF=1./rscc.pwfRef; % reverse prewhitening filter.
else
    dis.invPWF=1;  % inverse filter for particles
end;
mi.mbnOffset=rscc.membraneOffsetA/mi.pixA;

if ~isfield(mi,'ctf') || ~isfield(mi.ctf,'defocus')
    disp('No ctf information.');
    defoci=1;
else
    
    defoci=zeros(1,numel(mi.ctf));
    for i=1:numel(defoci)
        defoci(i)=mi.ctf(i).defocus;
    end;
end;
vesOk=all(mi.vesicle.ok,2);
if sum(vesOk)>0
    vesAmp=median(mi.vesicle.s(vesOk,1));
else
    vesAmp=0;
end;
% disp(['Defocus ' num2str(defoci,4) '  Dose ' num2str(mi.doses,3)...
%     '  Vesicle amp ',num2str(vesAmp*1000,3)]);

% screenSize=get(0,'screensize');
% Set the display size according to the preprocessing size--
% This avoids scaling bugs

dis.ds=mi.imageSize(1)/dis.ndis(1); % downsampling factor from original micrograph
if ~isfield(mi,'particle')
    mi.particle=struct;
end;

if dis.readAutopickPars && isfield(mi.particle,'autopickPars') && numel(mi.particle.autopickPars)>2 ...
        && mi.particle.autopickPars(1)>0
%     && ~dis.roboFitStep % Has valid parameters  % we'll allow this also
%     for robo-fitting.
    npars=numel(mi.particle.autopickPars);        
    if dis.readAutopickPars==1  % set the 'ap' parameters
            mask=false(1,max(npars,10));
            mask([1 2 10])=1;
            mask(npars+1:end)=[];
            dis.pars(mask)=mi.particle.autopickPars(mask);
            npars=sum(mask);
    else
            dis.pars(1:npars)=mi.particle.autopickPars;
    end;
    disp([num2str(npars) ' autopick parameters set.']);    
end;

% % warning('off');
% % save(dis.datName,'dis');  % store the preliminary settings for next time.
% save(dis.datName,'dis','-nocompression', '-v7.3');  % fast saving
% % warning('on');

% Set up the spectFactor value for this micrograph, depending on defocus
sFactor=1;
defocus=mi.ctf(1).defocus;
if dis.useSpectrumCorrectionTable % let it depend on defocus
    if defocus<dis.spectTable(1,1)
        sFactor=dis.spectTable(1,1);
    elseif defocus>dis.spectTable(end,1);
        sFactor=dis.spectTable(end,1);
    else
        sFactor=interp1(dis.spectTable(:,1),dis.spectTable(:,2),defocus);
    end;
    %     disp(['sFactor: ' num2str(sFactor)]);
end;
dis.pars(11)=sFactor;

ampFactor=1;
if dis.useAmpCorrectionTable % let it depend on defocus
    if defocus<dis.ampTable(1,1)
        ampFactor=dis.ampTable(1,1);
    elseif defocus>dis.ampTable(end,1);
        ampFactor=dis.ampTable(end,1);
    else
        ampFactor=interp1(dis.ampTable(:,1),dis.ampTable(:,2),defocus);
    end;
end;
dis.pars(12)=ampFactor;

% mag=dis.ndis(1)/size(rscc.mxCC,1);  % should be 1
m0=DownsampleGeneral(mc,dis.ndis);  % scale down merged image
% ---the following is code left over from not using the preprocessor image
% sizes for display:
% rscc.mxCC=DownsampleGeneral(rscc.mxCC,dis.ndis);
% rscc.mxVars=DownsampleGeneral(rscc.mxVars,dis.ndis);
% ne=size(rscc.eigenImgs,1);
% nei=NextNiceNumber(ne*mag);
% rscc.eigenImgs=DownsampleGeneral(rscc.eigenImgs,nei,mag,1);
% if numel(dis.invPWF)>1
%     dis.invPWF=DownsampleGeneral(dis.invPWF,nei,mag);
% end;
% rscc.mxVesInds=DownsampleNearest(rscc.mxVesInds,dis.ndis);
% rscc.mxTemplInds=DownsampleNearest(rscc.mxTemplInds,dis.ndis);
% rscc.mxRsos=DownsampleNearest(rscc.mxRsos,dis.ndis);
rs.blanks=rscc.mxVars==0;
% rs.blanks=DownsampleNearest(rscc.mxVars==0,dis.ndis);
rs.mCCU=rscc.mxCCU;
rs.mVesInds=rscc.mxVesInds<1;

if ~isfield(mi.vesicle,'ok')
    mi.vesicle.ok=[];
end;

% Now try to get the subtracted image
vesOk=false;
if dis.readVesicleImage
    vName=[mi.basePath 'Vesicles/' mi.baseFilename 'v.mrc'];
    %     [vName,vesOk]=CheckForImageOrZTiff(vName);
    [vName,vesOk]=CheckForAltImage(vName,sufExts);
    if vesOk
        mVesGood=DownsampleGeneral(ReadEMFile(vName),dis.ndis);
        mVesBad=0*mVesGood;
    end;
elseif dis.readSubImage
    [mv0,mvName,vesOk]=meReadMergedImage(mi,0,'vs');
%     mvName=[mi.basePath mi.procPath mi.baseFilename 'mv.mrc'];
    %     [mvName,vesOk]=CheckForImageOrZTiff(mvName);
%    [mvName,vesOk]=CheckForAltImage(mvName,sufExts);
    if vesOk
        disp(mvName);
        mv=DownsampleGeneral(mv0,dis.ndis);
        mVesGood=m0-mv;
        mVesBad=0*mVesGood;
    end;
end;
if ~vesOk % didn't already read it, check in rscc structure.
    if isfield(rscc,'mVesGood') && ~forceModelVesicles
        mVesGood=DownsampleGeneral(rscc.mVesGood,dis.ndis);
        mVesBad=DownsampleGeneral(rscc.mVesBad,dis.ndis);
    else
        % Create the vesicle models
        %         First, fix up the vesicle.ok array
        if size(mi.vesicle.ok,1)<numel(mi.vesicle.x)
            mi.vesicle.ok=true(numel(mi.vesicle.x),4);
        end;
        if size(mi.vesicle.ok,2)<4  % handle old mi files; all ves are good.
            mi.vesicle.ok=repmat(mi.vesicle.ok,1,4);
        end;
        
        if numel(mi.vesicle.ok)>0
            goodVes=all(mi.vesicle.ok(:,2:3),2);  % vesicle in range and refined
            badVes=mi.vesicle.ok(:,1) & ~goodVes; % refined but not in range
        else
            goodVes=[];
            badVes=[];
        end;
        disp(['Making Vesicles: ' num2str([sum(goodVes) sum(badVes)])]);
        
        mVesGood=meMakeModelVesicles(mi,dis.ndis,find(goodVes));
        mVesBad=meMakeModelVesicles(mi,dis.ndis,find(badVes));
    end;
end;

pixA=mi.pixA*mi.imageSize(1)/dis.ndis(1); % pixel size of displayed images

mVes=mVesGood+mVesBad;
rawMask=~meGetMask(mi,dis.ndis);

disp(['Defocus ' num2str(defoci,4) '  Dose ' num2str(mi.doses,3)...
    '  Vesicle amp ',num2str(vesAmp*1000,3)]);
disp([num2str(numel(mi.vesicle.x)) ' vesicles.']);


% Copy image and vesicle model into rscc
rscc.m0=m0;
rscc.mVes=mVes;

% Compute local variance
fc1=pixA/200;  % lowpass
fc2=pixA/1000; % baseline
fc3=pixA/300;  % var averaging
dev=GaussFilt(m0-mVes,fc1)-GaussFilt(m0-mVes,fc2);
if dis.rsMode
    rscc.mxVars=1000*sqrt(GaussFilt(dev.^2,fc3));
end;
rs.mxVars=DownsampleGeneral(rscc.mxVars,dis.ndis);

% All the images for display. 1 is unsub, 2 is vesicle-subtracted
imgs=zeros([dis.ndis 7],'single');
[imgs(:,:,1:2),dis.mulr,dis.addr]=rspFilterAndScaleImages(mi,dis,rscc);
imgs(:,:,3)=0;  % to receive model picked particles
imgs(:,:,4)=imgs(:,:,1); % to receive difference orig - model
imgs(:,:,5)=rscc.mxCCU*dis.ccuScale*100;
imgs(:,:,6)=imscale(rscc.mxTemplInds);
imgs(:,:,7)=imscale(rs.mxVars);
imgs(:,:,8)=dis.mulr*rscc.mVes+dis.addr;
imgs(:,:,9)=0;  % to receive overlap mask
dis.imgLabels{1}='';
dis.imgLabels{2}='Vesicles subtracted';
dis.imgLabels{3}='Model particles';
dis.imgLabels{4}='Image-model';
dis.imgLabels{5}='CC map';
dis.imgLabels{6}='Template indices';
dis.imgLabels{7}='Local variance';
dis.imgLabels{8}='Vesicle model';
dis.imgLabels{9}='Blank';

% masks are:  1: good vesicle ghosts;  2: bad vesicle ghosts;
%             3: variance above thresh 4: blank region
%             5: vesicle overlap
masks=zeros([dis.ndis 5],'single');
masks(:,:,1)=imscale(max(-mVesGood,0),1,1e-3);
if min(mVesBad(:))<0 % There are vesicles present
    masks(:,:,2)=imscale(max(-mVesBad,0),1,1e-3);
else
    masks(:,:,2)=0;
end;
if dis.rsMode
    masks(:,:,3)=((rs.mxVars>dis.pars(3)) | rawMask);
    masks(:,:,4)=rs.blanks;
else
    masks(:,:,3)=rawMask;
end;
dis.org=[0 0];
dis.mode=min(2,dis.mode);
dis.jpegCounter=1;
success=true;
end

function dis=ChangeTheWindowSize(dis,rscc)
rsz=size(rscc.mxCC);
if any(rsz)==0
    return
end;
if any(rsz~=dis.ndis)  % have to redraw the figure
    dis.ndis=rsz;
    dis.size=(dis.ndis);
    figure(1);
    if dis.clearFigure
        clf;
    end;
    % Keep the window in the old position
    pos=get(gcf,'position');
    set(gcf,'position',[pos(1:2) dis.size],'toolbar','none','resize','off');
    % Main display
    axsiz=dis.size(1)-3;
    aysiz=dis.size(2)-3;
    dis.ax1=axes('units','pixels','position',[2 3 axsiz aysiz]); %,'ticklength',[0 0]);
    dis.ax2=axes('position',[.8 0 .2 .2]);
    dis.ax3=axes('outerposition',[.8 0 .2 .2]);
    axes(dis.ax1);
end;
end