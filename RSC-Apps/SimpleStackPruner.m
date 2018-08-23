% SimpleStackPruner.m
%  Based on SimpleRSPicker, show picked particles in context and allow
%  particles to be removed from a stack.

% skipToLoop=exist('si','var') && exist('stack','var') && exist('dis','var');
skipToLoop=0;

if ~skipToLoop
version=12;
extraDs=1;  % 2 for retina display
showDualStacks=1;
% Retrieve parameters from a file in the program directory
pa=fileparts(which('SimpleStackPruner'));
datName=[AddSlash(pa) 'SimpleStackPrunerDat.mat'];
disOk=0;
if exist(datName,'file')
    dis=load(datName);
    dis=dis.dis;
    dis.datName=datName;
    if isfield(dis,'version')
        disOk=dis.version==version && dis.siValid;
        disp('Loading settings');
    end;
end;
newDis=0;
if ~disOk
    disp('Loading defaults');
    dis=struct;
    dis.version=version;
    dis.datName=datName;
    dis.stackds=1; % downsampling for stack display.  4 is standard
    dis.nCrop=48;
    dis.showPW=0;
    dis.readVesicleImage=true; % attempt to read vesicle image from file
    dis.imageMode=false;
    
    % Box colors
    bColorErased=[.3 .3 .3];
    bColorMan=[.7 .7 .2];  % yellow
    bColorManNotRefined=[.7 .6 .3]; %orangish
    bColorAuto=[.4 .7 .2];  % greenish
    bColorBkd=[1 .7 .2]; % bright orange-- picked particle
    bColorVesicle=[.2 0 1]; % blue vesicle
    bColorBadVesicle=[.5 0 0]; % red bad vesicle
    dis.boxColors=[bColorErased; bColorMan; bColorAuto; bColorBkd;
        bColorVesicle; bColorManNotRefined; bColorBadVesicle];
    dis.lineWidth=2;
    dis.corners=[1  1  1 1 .41 1.2 .41
        0 .8 .8 1 .41 1.2 .41]'; %  'eroded' box corners
    
    dis.markers=['ys';'cs';'gs';'ms'];
    dis.clickMarker=['gs';'ks'];
    
    dis.showDualStacks=0;
    dis.stackStride=1;
    
    % overlay colors
    dis.maskColor=[.9 .75 .8];
    dis.blankColor=[1 .85 .8];
    dis.overlapColor=[.95 .8 1];
    dis.ghostColor=[.7 .7 1];
    dis.ghostColorBad=[1 .7 .6];
    dis.ghostAmp=.6;
    
    % initial display
    dis.polarity=1;
    dis.ndis=[768 768]; %%
    dis.size=[768 768];
    dis.size=min(dis.size,dis.ndis); % Force the display to be no bigger than the image
    dis.org=[0 0]; %pixels run from org+1 to org+size
    
    %     default parameters
    dis.mode=1;  % vesicles not subtracted
    dis.stackStep=1;
    dis.umul=.6;
    dis.showMask=0;
    dis.showGhosts=0;
    dis.showBoxes=1;
    dis.listParticleInfo=1;
    dis.contrast=[1 1 .001];
    dis.varThresh=40;
    dis.pars=[.55 1 40 35 100 120 120 200];  % min, max amp; max var; rso offset;
    %     particle blank radius, vesicle blank radius, maxBob, border.
    dis.pars(20)=200;  % box size.
    dis.minDist=dis.pars(20);  % distance in original pixels, based on box size
    dis.filter=[1000 20 0];  % inverse frequency in A
    
    dis.basePath='./';
    dis.siPath='';
    dis.infoName='';
    dis.stackIndex=1;
    dis.finished=0;
    dis.siValid=0;
    newDis=1;
end;
%%
% Try to load the previous si file

if exist(dis.siPath,'dir')
    cd(dis.siPath);
end;

% Set the first command
if (dis.finished || newDis || ~dis.siValid) % We need to open a new file
    b='o';  % ask for a new file
else
    b='v';  % 'revert' to latest file.
end;
oldB=0;

% % disp('Make figure');
% screenSize=get(0,'screensize');
% xsiz=min(screenSize(3),dis.size(1)+3);
% ysiz=min(screenSize(4)-50,dis.size(2)+3);
% figure(1);
% clf;
% % Put the window near the top middle of the screen
% set(gcf,'position',[(screenSize(3)-dis.size(1))/2 ((screenSize(4)-50)-ysiz)*0.9 xsiz ysiz],'toolbar','none','resize','off');
% % Main display
% axsiz=xsiz-3;
% aysiz=ysiz-3;
% ax1=axes('units','pixels','position',[2 3 axsiz aysiz]); %,'ticklength',[0 0]);
% ax2=axes('position',[.8 0 .2 .2]);
% ax3=axes('outerposition',[.8 0 .2 .2]);
% axes(ax1);

refreshReconstruct=0;  % flag to update the reconstruction display
stackIndex=1;
mIndexOld=0;
dis.currentBoxSize=dis.pars(20);
prevPartIndex=0;

else
    disp('Restarting...');
    b=0;
    disp('Init display');
    nim=InitDisplay(dis,stack,markInds);
    disp([num2str(nim) ' particles.']);
end;  % ~ skip to loop
%%
% % interactive loop
while (b~='q') && (b~='Q') && (b~=4) % q = quit; double-click=quit.
    
    switch b
        
        case {'o'}  % open a stack file

            [dis,si,stack,markInds]=rspLoadActiveStacks(dis);
            n=size(stack,1);
            if n<1
                disp('No stack loaded.');
            end;
%             if showDualStacks
%                 [si0, stack0, names, siPath, ustack0]=reLoadStackFiles;
%                 dis.stackStep=2;
%             else
%                 [si0, stack0, names, siPath]=reLoadStackFiles;
%                 dis.stackStep=1;
%             end;
%             disp('loaded:');
%             disp(names);
%             j=size(si0.activeFlags,2);
%             if (b=='o') && j>1 % possibly trim the activeFlags
%                 afIndex=PickActiveFlagSet(si0);
% %                 disp([num2str(size(si0.miIndex,1)) ' particles total.']);
% %                 for i=1:j
% %                     if numel(si0.activeFlagLog)<i
% %                         si0.activeFlagLog{i,1}='';
% %                     end;
% %                     disp([num2str([i sum(si0.activeFlags(:,i))]) '  ' si0.activeFlagLog{i}]);
% %                 end;
% %                 j=MyInput('Start with data at selection ',j);
% %                 if j<1 % reset the active flags
% %                     si0.activeFlags=true(size(si0.miIndex));
% %                     j=1;
% %                     disp('Resetting the active flags.');
% %                 end;
% %                 afIndex=j;
%                 active=si0.activeFlags(:,afIndex);
%             else
%                 active=si0.activeFlags(:,1);
%                 afIndex=1;
%             end;
%             q=input('Show only these particles? ','s');
%             if lower(q)=='y'
%               if showDualStacks
%                   [si,stack,ustack]=rsStackSplit(si0.activeFlags(:,afIndex),si0,stack0,ustack0);
%               else
%                   [si,stack]=rsStackSplit(si0.activeFlags(:,afIndex),si0,stack0);
%               end;
%               afIndex=PickActiveFlagSet(si);
%               active=si.activeFlags(:,afIndex);
%             else
%                 si=si0;
%                 stack=stack0;
%                 if showDualStacks
%                     ustack=ustack0;
%                 end;
%             end;
%             nim=size(stack,3);
%             if showDualStacks % we'll interleave them
%                 stack0=stack;
%                 usd=std(ustack(:));
%                 stack(:,:,1:2:2*nim)=dis.umul*(ustack-1.5*usd);
%                 stack(:,:,2:2:2*nim)=stack0;
%                 dis.stackStep=2;
%             else
%                 dis.stackStep=1;
%             end;
%             n=size(stack,1);
            nim=InitDisplay(dis,stack,markInds);
            active=false(nim,1);
            for i=1:numel(markInds)
                active(markInds{i})=true;
            end;
            
            
%             size(stack,3)/dis.stackStride;
%             disp([nim ' images.']);
%  
%             figure(2);
%             imagsar(dis.polarity*BinImage(Crop(stack,dis.nCrop,1),dis.stackds/extraDs),dis.contrast(3),0);
% 
%             disp('marking particles');
%             for k=1:numel(markInds)
%                 mk.marker=dis.markers(k,:);
%                 mk.markerSize=n/dis.stackds;
%                 imagsar('mark',markInds{k},mk);
%             end;

        case {1 2 3 'x'}  % mouse button  1: select; 3: mark bad; 2: flip range bad
%             key 'x': mark the whole range bad
%               you select a range by left-clicking the the first element,
%               then button 2 or x on the last element.
            dis.stackIndex=ceil(clickIndex/dis.stackStride);
            disp([' particle ' num2str(dis.stackIndex)]);
            %             stackIndex=MyInput(' Particle index',stackIndex);
            if dis.stackIndex>0 && dis.stackIndex <=numel(si.miIndex)
                if dis.imageMode
                    figure(1);
                    mIndex=si.miIndex(dis.stackIndex);
                    pIndex=si.miParticle(dis.stackIndex);
                    
                    % load the image
                    mi=si.mi{mIndex};
                    mi.basePath='';
                    if mIndex~=mIndexOld
                        disp(['Loading ' mi.baseFilename]);
                        mc=meReadMergedImage(mi);
                        if dis.readVesicleImage
                            vName=['Vesicles/' mi.baseFilename 'v.mrc'];
                            vOk=exist(vName,'file');
                            if vOk
                                mv=ReadMRC(vName);
                                rscc.mVes=DownsampleGeneral(mv,dis.ndis);
                            end;
                        end;
                        if numel(mc)<2  % unsuccessful read
                            mc=single(zeros(1024));
                            disp('No merged image found.');
                        elseif dis.showPW
                            fHP=.005;
                            nZeros=1;
                            ds=si.pixA/mi.pixA;
                            nb=size(si.ctfs,1);
                            H=meGetNoiseWhiteningFilter(mi,nb,ds,nZeros,fHP);  % for micrograph
                            %                        Get all the particles in this image and compute
                            %                        the spectrum
                            partInds=find(si.miIndex==mIndex);
                            partImgs=stack(:,:,partInds);
                            pmsk=1-fuzzymask(nb,2,nb*0.4,nb*.04);
                            S=RimSpectrum(partImgs,pmsk);
                            medS=median(S(nb/4:3*nb/8));
                            S=S/medS;
                            %                        S(round(nb*0.4):end)=NaN;
                            figure(3);
                            fs=(0:nb/2-1)/(nb*si.pixA);
                            plot(fs,[sectr(H) sqrt(S)]);
                            xlabel('Frequency, A^{-1}');
                            ylabel('PW Filter');
                        end;
                        
                        dis.ds=mi.imageSize(1)/dis.ndis(1);
                        
                        m0=DownsampleGeneral(mc,dis.ndis);
                        rscc.m0=m0;
                        imgs=zeros([dis.ndis 7],'single');  % Images for display
                        imgs(:,:,1:2)=rspFilterAndScaleImages(mi,dis,rscc);
                        masks=zeros([dis.ndis,5],'single');
                        dis.org=[0 0];
                        dis.mode=1;
                        mIndexOld=mIndex;
                        oldPIndex=0;
                    end;
                    if oldPIndex>0
                        mi.particle.picks(oldPIndex,3)=oldPFlag;
                    end;
                    oldPFlag=mi.particle.picks(pIndex,3);
                    mi.particle.picks(pIndex,3)=48;  % mark as background
                    oldPIndex=pIndex;
                    
                    [picks, ptrs]=rspLoadPicksFromMi(mi);
                    %                 [picks, ptrs]=rspDeleteBadAutoPicks(dis,picks,ptrs);
                    refreshReconstruct=1;
                    if ~any(rspCountParticles(picks))
                        disp('No particles found');
                    end;
                    %                 disp([num2str(rspCountParticles(picks)) ' particles loaded']);
                    figure(1);
                    rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                    set(gcf,'name',mi.baseFilename);
                    figure(2);
                end;  % if imageMode
 % -----------handle mouse buttons ---------------               
                if b==1  % toggle active
                    if active(dis.stackIndex)                    
                       mk.marker=dis.clickMarker(2,:); % erase
                       active(dis.stackIndex)=false;
                    disp('not active.');
                    else
                        mk.marker=dis.clickMarker(1,:); % mark active
                        active(dis.stackIndex)=true;
                    disp('active.');
                    end;
                    mk.markerSize=n/dis.stackds;
                    imagsar('mark',dis.stackIndex*dis.stackStride,mk);
                    prevPartIndex=dis.stackIndex;
                elseif (b==2 || b=='x') && prevPartIndex>0
                    q=active(prevPartIndex); % shift-click: extend state;
                    for p=dis.stackStride*(prevPartIndex+1:dis.stackIndex)
                        if q
                            mk.marker=dis.clickMarker(1,:); % activate
                            active(p)=true;
                        else
                            mk.marker=dis.clickMarker(2,:); % erase
                            active(p)=false;
                        end;
                        mk.markerSize=n/dis.stackds;
                        imagsar('mark',p,mk);
                    end;
                elseif b==3
                    if active(dis.stackIndex)
                        mk.marker='ro';
                    else
                        mk.marker='ko';
                    end;
                    mk.markerSize=0.6*n/dis.stackds;
                    imagsar('mark',clickIndex,mk);
                    active(dis.stackIndex)=~active(dis.stackIndex);
                end;
                
                %                 if b==2
                %                     if active
                %                     active(partIndex)=~active(partIndex);
                %                     markedBad(partIndex)=~markedBad(partIndex);
                %                 end;
                %                 if active(partIndex)
                %                     if markedBad(partIndex)
                %                         mk.marker='ko';
                %                         mk.markerSize=n*.9/dis.stackds;
                %                         markedBad(partIndex)=false;
                %                     else
                %                     mk.marker='bs';
                %                     mk.markerSize=n/dis.stackds;
                %                     end;
                %                 else
                %                     mk.marker='ro';
                %                     mk.markerSize=n*.9/dis.stackds;
                %                 end;
            end; % if stackIndex
        case 'B'
            dis.pars(20)=MyInput(' Box size, A ',dis.pars(20));
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'c'  % set contrast
            disp('Setting contrast');
            if dis.imageMode
                dis.contrast(1)=MyInput('Black contrast',dis.contrast(1));
                dis.contrast(2)=MyInput('White contrast',dis.contrast(2));
                imgs(:,:,1:2)=rspFilterAndScaleImages(mi,dis,rscc);
                figure(1);
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            end;
            figure(2);
            dis.contrast(3)=MyInput('Imagic contrast',dis.contrast(3));
        case 'd'  % display mode (toggle 1-2, 2-3, 1-end)
            if dis.imageMode
                dis.mode=dis.mode+1;
                if dis.mode>2 % lower case toggles norm / ves-subtr
                    dis.mode=1;
                end;
                if dis.mode==2 && sum(sum(imgs(:,:,2)))==0  % no vesicle model
                    disp(' making vesicles...');
                    if size(mi.vesicle.ok,1)<numel(mi.vesicle.x)
                        mi.vesicle.ok=true(numel(mi.vesicle.x),4);
                    end;
                    mVes=meMakeModelVesicles(mi,dis.ndis,find(mi.vesicle.ok(:,1)));  % Actual vesicle model
                    rscc.mVes=mVes;
                    imgs(:,:,1:2)=rspFilterAndScaleImages(mi,dis,rscc);
                    disp(' done.');
                    %             elseif b=='e' && dis.mode>3
                    %                 dis.mode=2;
                    %             end;
                    %             if dis.mode > size(imgs,3) % upper case cycle through all imgs
                    %                 dis.mode=1;
                    %             end;
                    %             if dis.mode==3 && refreshReconstruct
                    %                 imgs(:,:,3)=rspReconstructParticles(dis,mi,picks,ptrs,rscc);
                    %                 refreshReconstruct=0;
                    %             end;
                end;
                figure(1);
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                figure(2);
            end;
            q
        case 'f'  % set filtering
            if dis.imageMode
                disp('Setting the filter:');
                dis.filter(1)=MyInput('Highpass filter, A',dis.filter(1));
                dis.filter(2)=MyInput('Lowpass filter, A',dis.filter(2));
                imgs(:,:,1:2)=rspFilterAndScaleImages(mi,dis,m0,mVes);
                figure(1);
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                figure(2);
            end;
        case 'P' % set major params
            dis.imageMode=MyInput('Image mode',dis.imageMode);
            
    end; % switch
    oldB=b;  % store the previous key
    pause(.05);
    [clickIndex,coords,b]=imagsar;
    if numel(b)<1
        b=0;
    end;
%      if b, disp(b); end;
    %     if b~=0
    %         b
    %     end;
    %     axes(ax1);
end;  % while b~='q'
%
originalNo=numel(active)
finalNo=sum(active)
si1=si;

si1.activeFlags(:,end+1)=active;
si1.activeFlagLog{end+1,1}=[date '  SimpleStackPruner'];
% si1=rsStackSetActiveFlags(active,si0,'SimpleStackPruner');

str=input('Write the si file y: add new flags; r: replace flags; n: don''t write [n]? ','s');
if numel(str)>0 && lower(str(1))=='y'
    si1.activeFlags(:,end+1)=active;
    si1.activeFlagLog{end+1,1}=[date '  SimpleStackPruner'];
    siTemp=si;
    si=si1;
    name=names{end};
    name=[siPath name];
    save(name,'si');
    disp([name ' written.']);
elseif numel(str)>0 && lower(str(1))=='r'
    si1.activeFlags(:,afIndex)=active;
    si1.activeFlagLog{afIndex,1}=[date '  SimpleStackPruner'];
    siTemp=si;
    si=si1;
    name=names{end};
    name=[siPath name];
    save(name,'si');
    disp([name ' written.']);
else
    disp('Nothing written.');
end;

% save(datName,'dis');
% b=0;


function nim=InitDisplay(dis,stack,markInds)
            nim=size(stack,3)/dis.stackStride;
            n=size(stack,1);
            nd=dis.nCrop/dis.stackds;
            figure(2);
            disp('downsampling');
            msk=Gaussian(nd,2,nd*.2);
            if dis.stackds>1
                dsStack=Downsample(Crop(stack,dis.nCrop,1),nd,1,msk);
            else
                dsStack=GaussFilt(stack,.2,1);
            end;
            disp('display');
            imagsar(dis.polarity*dsStack,dis.contrast(3),0);

            disp('marking particles...');
            n=size(stack,1);
            for k=1:numel(markInds)
                k1=min(k,size(dis.markers,1));
                mk.marker=dis.markers(k1,:);
                mk.markerSize=n/dis.stackds;
                imagsar('mark',dis.stackStride*markInds{k},mk);
            end;
            disp('done.');
end