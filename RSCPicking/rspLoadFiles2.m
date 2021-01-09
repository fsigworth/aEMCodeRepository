function [dis, mi, rscc, rs, imgs, masks, rawMask,success]=rspLoadFiles2(dis,si)
% Given the info file name and path (in the dis structure), load the mi and rscc
% files.
% 16-Dec-20: This is version 2 where I have eliminated code that doesn't
% assume the existence of a modern rscc file and loads images from elsewhere.
% I also include offsets for mi.useMicrographCoords. We add the cropOffset
% here, and remove it in rspStorePicksInMi.

% Default return values
mi=[];
rscc=[];
rs=[];
imgs=[];
masks=[];
rawMask=[];
success=0;

dis.basePath=AddSlash(pwd);
disp(' ');
if dis.miIndex
    disp(['file index: ' num2str(dis.miIndex)]);
end;

dis.miValid=0;
if dis.miMode
%     miName=[AddSlash(dis.infoPath) dis.infoName];
    miName=[dis.infoName];
    if exist(miName,'file')
        mi=ReadMiFile(miName,1);  % get the mi structure
        disp([miName ' loaded.']);
        dis.miValid=1;
    else
        disp(['mi file not found: ' miName]);
        return
    end;
else
    mi=si.mi{dis.miIndex};
    dis.miValid=1;
end;
mi.basePath=dis.basePath;
set(gcf,'name',[dis.infoName ' - loading...']);

% Deal with original micrograph coordinates, if we're using them.
% This will allow us to do particle extraction using Relion.
if ~isfield(mi,'padImageSize')
    mi.padImageSize=NextNiceNumber(mi.imageSize,5,8);
end;
if dis.forceMicrographCoords
    mi.useMicrographCoords=1;
end;
if isfield(mi,'useMicrographCoords') && mi.useMicrographCoords
    %     We'll modify the particle and vesicle coordinates to account for the
    %     offset from padding
    cropOffset=floor((mi.padImageSize-mi.imageSize)/2);
    if isfield(mi,'particle') && isfield(mi.particle,'picks') ...
            && size(mi.particle.picks,1)>0 && size(mi.particle.picks,2)>0
        %         we offset all the particle coordinates
        mi.particle.picks(:,1)=mi.particle.picks(:,1)+cropOffset(1);
        mi.particle.picks(:,2)=mi.particle.picks(:,2)+cropOffset(2);
    end;
    %     Offset the vesicle coordinates and stash the originals.
    dis.origVesicle.x=mi.vesicle.x;
    dis.origVesicle.y=mi.vesicle.y;
    mi.vesicle.x=mi.vesicle.x+cropOffset(1);
    mi.vesicle.y=mi.vesicle.y+cropOffset(2);
end;

% handle bug where weights are zero.
if sum(mi.weights)==0
    mi.weights=ones(numel(mi.ctf),1);
end;

rawMask=~meGetMask(mi,dis.ndis);

% ---------Pick up the rscc file and images---------
if isfield(mi,'procPath_sm') % Look in this folder
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
    nRscc=size(rscc.mxCC);
    if any(nRscc ~= dis.ndis)
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
    return
end;

% Now that display size has been established to match rscc images...
dis.ds=mi.padImageSize(1)/dis.ndis(1); % downsampling factor from original micrograph
dis.pixA=mi.pixA*dis.ds;

% pick up the raw and subtracted images from the rscc file.
if ~isfield(rscc,'m0') % We can't read the images directly
    rscc=LoadImageFiles(mi,rscc); % handle old-form files
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
            
        defoci=zeros(1,numel(mi.ctf));
        for i=1:numel(defoci)
            defoci(i)=mi.ctf(i).defocus;
        end;
    vesOk=all(mi.vesicle.ok,2);
    if sum(vesOk)>0
        vesAmp=median(mi.vesicle.s(vesOk,1));
    else
        vesAmp=0;
    end;
    
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
    
    rs.blanks=rscc.mxVars==0;
    % rs.blanks=DownsampleNearest(rscc.mxVars==0,dis.ndis);
    rs.mCCU=rscc.mxCCU;
    rs.mVesInds=rscc.mxVesInds<1;
    
%     if ~isfield(mi.vesicle,'ok')
%         mi.vesicle.ok=[];
%     end;
%     
    
    rawMask=~meGetMask(mi,dis.ndis);
    
    disp(['Defocus ' num2str(defoci,4) '  Dose ' num2str(mi.doses,3)...
        '  Vesicle amp ',num2str(vesAmp*1000,3)]);
    disp([num2str(numel(mi.vesicle.x)) ' vesicles.']);
    
    
    % Copy image and vesicle model into rscc
    rscc.mVes=rscc.m0-rscc.m1;
    mVesGood=rscc.mVes;
    mVesBad=0*mVesGood;
    
    % Compute local variance
    fc1=dis.pixA/200;  % lowpass
    fc2=dis.pixA/1000; % baseline
    fc3=dis.pixA/300;  % var averaging
    dev=GaussFilt(rscc.m0-rscc.mVes,fc1)-GaussFilt(rscc.m0-rscc.mVes,fc2);
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
    dis.imgLabels= ...
        {'';
        'Vesicles subtracted';
        'Model particles';
        'Image-model';
        'CC map';
        'Template indices';
        'Local variance';
        'Vesicle model';
        'Blank'};

%     
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
        masks(:,:,3)=((rs.mxVars>dis.pars(3)) | rawMask);
        masks(:,:,4)=rs.blanks;
    dis.org=[0 0];
    dis.mode=min(2,dis.mode);
    dis.jpegCounter=1;
    success=true;
end % -------------------Main----------------





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

    
    
    function rscc=LoadImageFiles(mi,rscc)
%     Code for getting image and subtracted image for old-style rscc files
        %sufExts={'s.mrc' 'z.tif' '.mrc'}; % suffix and extension options for small, compressed or full files.
    % drawnow;
    [rscc.m0, mergedName,ok]=meReadMergedImage(mi,0,'s');
    %[mergedName,ok]=CheckForAltImage([mi.procPath mi.baseFilename 'm.mrc'],sufExts);
    % [mc,mergedName]=meReadMergedImage(mi);
    % mc=ReadEMFile(mergedName);
    disp(mergedName);
    if ~ok  % unsuccessful read
        rscc.m0=single(zeros(960));
        disp('No merged image found.');
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
    
    pixA=mi.pixA*dis.ds; % pixel size of displayed images
    
    mVes=mVesGood+mVesBad;
    rawMask=~meGetMask(mi,dis.ndis);
    
    disp(['Defocus ' num2str(defoci,4) '  Dose ' num2str(mi.doses,3)...
        '  Vesicle amp ',num2str(vesAmp*1000,3)]);
    disp([num2str(numel(mi.vesicle.x)) ' vesicles.']);
    
    
    % Copy image and vesicle model into rscc
    rscc.m0=m0;
    rscc.mVes=mVes;
    
    end
