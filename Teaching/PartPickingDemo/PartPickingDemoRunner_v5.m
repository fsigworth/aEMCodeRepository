% PartPickingDemoRunner_5.m

nonScale=1.2;
% pixA=.822;

load('TemplatePickerDat.mat'); % Get the po struct.
po.starInPath='CtfFind/job003/';  %%%%%%

blankRadiusA=[100 100]; % Blank for non-particles, particles
        msz=30; % marker size

th0=.4;  % K2 starting value
% th0=100; % K3 starting value
displayFilter=.2;
maxThresholdMul=2;
interpolateImage=2;

% Read the star file
% cd /Users/fred/EMWork/Yangyu/20190325/MotionCorr/job061_micrographs_ctf_Topaz
% cd  /Volumes/Drobo4/Yangyu/20200603/
outPath='TemplatePicker4_stars/';
matInPath=po.matOutPath;

CheckAndMakeDir(outPath,1);

[~,d]=ReadStarFile([po.starInPath 'micrographs_ctf.star']);
d=d{1};
nl=numel(d.rlnMicrographName);

ind=1;

% Get the first mat file
[pa,baseName,ex]=fileparts(d.rlnMicrographName{ind});
matName=[matInPath baseName '.mat'];
load(matName);

ds=po.ds;
blankRadius=blankRadiusA/(pixA*ds); % bad, good masks

blankMasks=cell(2,1);
blankMasks{1}=(1-fuzzymask(2*ceil(blankRadius(1)+3),2,blankRadius(1),2));
blankMasks{2}=(1-fuzzymask(2*ceil(blankRadius(2)+3),2,blankRadius(2),2));

CheckAndMakeDir(outPath,1);

figure(4);
MyBusyread('init');

b='N';
ind=1;
displayMode=1;
boxesOn=2;
interpFactor=interpolateImage;

while b~='q'
    doLoad=0;
    doSearch=1; % by default we search and update display
    switch b
        case 'u' % up the bad scaling
            nonScale=nonScale*1.05;
        case 'i' % down with bad scaling
            nonScale=nonScale/1.05;
        case 'j' % increase threshold
            th0=th0*1.05;
        case 'k' % keep more
            th0=th0/1.05;
            
        case {'n' ' ' 'R'} % Get next file or repeat
            ind=ind+1;
            if ind>nl
                beep;
                ind=nl;
            end;
            doLoad=1;
        case {'N' 'p'} % Get the previous file
            ind=ind-1;
            if ind<1
                beep;
                ind=1;
            end;
            doLoad=1;
        case 'g' % go to file index
            ind=MyInput('go to file index ',ind);
            ind=max(1,ind);
            doLoad=1;
        case 'b' % change box display
            boxesOn=boxesOn+1;
            if boxesOn>2
                boxesOn=0;
            end;
        case {'e' 'r'} % show expectation or residual
            displayMode=3-displayMode;  % either 1 or 2
        case 'z'    % zero out the particles
            picks=zeros(0,4);
            badPicks=zeros(0,4);
        otherwise
            doSearch=0;
    end;
    
    if doLoad
        if ind>nl
            disp('End of files');
            return
        end;
        [pa,baseName,ex]=fileparts(d.rlnMicrographName{ind});
        matName=[matInPath baseName '.mat'];
        if exist(matName,'file')
            load(matName);
            nGroups=size(po.ccImgs,3);
            mdmul=0; % display scaling not set
            displayMode=1;
        else
            beep;
            disp(' ');
            disp(['----The file ' matName ' wasn''t found.']);
        end;
    end;
    
%     Image display
    
    
    if doSearch
        if displayMode==1 % Always 1 if we haven't displayed yet
            
            % Set up scaling of the filtered image
            mdf=GaussFilt(po.md,displayFilter);
            [mds,mdmul,mdadd]=imscale(mdf,256,.001);
            imaga(Downsample(mds,interpolateImage*po.nd));
            interpFactor=interpolateImage;
        else % display mode 2: show residual or expectation
            np=size(picks,1);
            modelImage=zeros(nd,'single');
            for i=1:np
                modelImage=modelImage+ExtractImage(picks(i,3)*po.cTemplates(:,:,picks(i,4)),...
                    picks(i,1:2),nd,1);
            end;
            if b=='r'
                doSearch=0;
                modelImage=mdf-modelImage;
            end;
            %                  imaga(Downsample(mdmul*modelImage+mdadd,interpolateImage*po.nd));
            imaga(mdmul*modelImage+mdadd);
            interpFactor=1;
        end;
        
        
        % set up to search
        stackScales=ones(nGroups,1);
        stackScales(end)=nonScale;
        
        
        ccTotImgs=po.ccImgs(:,:,1);
        nd=size(ccTotImgs);
        ccTotInds=po.ccInds(:,:,1);
        ccGroups=ones(nd,'uint8'); % goo particles have group=1
        for i=2:nGroups
            ccXImgs=stackScales(i)*po.ccImgs(:,:,i);
            next=ccXImgs>ccTotImgs;
            ccTotImgs(next)=ccXImgs(next);
            nextInds=po.ccInds(:,:,i);
            ccTotInds(next)=nextInds(next);
            ccGroups(next)=i;  % bad particles will have group>1
        end;
        ccGoods=any(ccGroups==po.goods,3);
                
        ccMax=inf;
        
        ccImgWk=ccTotImgs;
        
        picks=zeros(0,4);
        badPicks=zeros(0,4);
        
        th1=5*th0;
        
        if b~='z' % We always search unless we've zeroed the picks.
            
            title('---searching-----');
            drawnow;
            
            hold on;
            while ccMax>th0
                [ccMax,i,j]=max2d(ccImgWk);
                if ccMax<th0
                    break
                end;
                good=ccGoods(i,j);
                %             mysubplot(221);
                %             imags(ccImgWk);
                %             mysubplot(222);
                %             imags(ccGoods);
                ccImgWk=Mask(ccImgWk,[i j],blankMasks{good+1});
                %             mysubplot(223);
                if ccMax<th1 % lies below the bad particle threshold
                    if good % a true particle
                        marker='gs';
                        plot(i*interpFactor,j*interpFactor,marker,'markersize',msz);
                        picks(end+1,:)=[i j double(ccMax) double(ccTotInds(i,j))];
                    elseif boxesOn>1
                        marker='ro';
                        plot(i*interpFactor,j*interpFactor,marker,'markersize',4);
                        badPicks(end+1,:)=[i j double(ccMax) double(ccTotInds(i,j))];
                    end;
                end;
                
            end; % while
            hold off;
        end; % if ~z
        
        txt=num2str([th0 nonScale size(picks,1)],4);
        text(1,1,txt,'verticalalignment','bottom','horizontalalignment','left',...
             'fontsize',16,'color',[1 1 0]);
        title(matName,'interpreter','none');
        drawnow;
        %     axis off;
        updateDisplay=0;
    end; % if doSearch
    
    % If we're repeating, any key stops.
    oldB=b;
    pause(0.2);
    [ix,iy,b]=MyBusyread;
    if numel(b)<1 || ~ischar(b)
        b=char(0);
    end
    if oldB=='R' && ischar(b) && b==char(0)
        b=oldB;
    elseif oldB=='R'
        disp('Robopicking stopped.');
    end;
    
    np=size(picks,1);
    
    if any(b==['nNqzR'])
        
        % Construct the output star file
        pk=struct;
        pk.rlnCoordinateX=picks(:,1)*ds-po.offsets(1)-1;
        pk.rlnCoordinateY=picks(:,2)*ds-po.offsets(2)-1;
        pk.rlnAutopickFigureOfMerit=picks(:,3);
        pk.rlnClassNumber=zeros(np,1);
        pk.rlnAnglePsi=zeros(np,1);
        outStarName=[outPath po.baseName '_autopick.star'];
        WriteStarFileStruct(pk,'',outStarName);
        drawnow;
        disp(['Index ' num2str(ind) ': Written: ' outStarName ' ' num2str(size(picks,1)) ' picks']);
        
        if b=='q' % For the last one, read back the coordinates and show the match.
            % Read back the coordiantes
            [nms,pk0]=ReadStarFile([outPath po.baseName '_autopick.star'],1);
            pk0=pk0{1};
            offsets=po.offsets;
            
            imags(mdf);
            hold on;
            % plot(pk.rlnCoordinateX*ds+offsets(1)+1,rlnCoordinateY*ds+offsets(2)+1),'gs');
            plot( (pk.rlnCoordinateX+offsets(1))/ds+1,(pk.rlnCoordinateY+offsets(2))/ds+1,'gs','markersize',msz);
            plot(badPicks(:,1),badPicks(:,2),'rs','markersize',4);
            hold off;
        end;
        
        
    end;
end; % while
disp('done');


return

%%

% Code to make figure for class

pk=picks(51,:);

templ=Crop(po.cTemplates(:,:,pk(4)),nd);
cc=real(ifftn(fftn(po.md).*conj(fftn(ifftshift(templ)))));
imaga(cc*300+100); hold on; plot(663,425,'gs','markersize',40); hold off;


%%

% Demonstrate signal in noise

N=randn(nd);
Nf=GaussFilt(randn(nd),.2);
si=Crop(po.cTemplates(:,:,17),nd);
X=Nf+5*si;

imags(X);

cc=real(ifftn(fftn(X).*conj(fftn(ifftshift(si)))));
% imags(cc);

%%
% Compare picks
[nms,pk0]=ReadStarFile([po.baseName '_autopick.star'],1);
pk0=pk0{1};

% imags(mdf);
hold on;
% plot(pk.rlnCoordinateX*ds+offsets(1)+1,rlnCoordinateY*ds+offsets(2)+1),'gs');
plot( (pk.rlnCoordinateX+offsets(1))/ds+1,(pk.rlnCoordinateY+offsets(2))/ds+1,'gs');
plot((pk0.rlnCoordinateX+offsets(1))/ds+1,(pk0.rlnCoordinateY+offsets(2))/ds+1,'rs');
hold off;
