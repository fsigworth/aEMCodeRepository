% PartPickingDemoRunner.m

nonScale=1.2;
pixA=.822;
ds=4;
blankRadius=[200 100]/(pixA*ds); % bad, good masks
th0=.4;
displayFilter=.2;

% Read the star file
cd /Users/fred/EMWork/Yangyu/20190325/MotionCorr/job061_micrographs_ctf_Topaz
outPath='OutStarFiles/';

[~,d]=ReadStarFile('micrographs_ctf.star');
d=d{1};
nl=numel(d.rlnMicrographName);

ind=1;

% Get the first mat file
[pa,baseName,ex]=fileparts(d.rlnMicrographName{ind});
matName=[baseName '.mat'];
load(matName);

blankMasks=cell(2,1);
blankMasks{1}=(1-fuzzymask(2*ceil(blankRadius(1)+3),2,blankRadius(1),2));
blankMasks{2}=(1-fuzzymask(2*ceil(blankRadius(2)+3),2,blankRadius(2),2));

CheckAndMakeDir(outPath,1);

figure(4);
MyBusyread('init');

b='N';
ind=1;
displayMode=1;
boxesOn=1;

while b~='q'
    doLoad=0;
    doSearch=0;
    switch b
        case 'u' % up the bad scaling
            nonScale=nonScale*1.05;
        case 'i' % down with bad scaling
            nonScale=nonScale/1.05;
        case 'j' % increase threshold
            th0=th0*1.05;
        case 'k' % keep more
            th0=th0/1.05;
        case {'n' ' '}
            ind=ind+1;
            if ind>nl
                beep;
                ind=nl;
            end;
            doLoad=1;
        case 'N'
            ind=ind-1;
            if ind<1
                beep;
                ind=1;
            end;
            doLoad=1;
        case 'b'
            boxesOn=boxesOn+1;
            if boxesOn>2
                boxesOn=0;
            end;
        case {'e' 'r'}
            displayMode=3-displayMode;  % either 1 or 2
        case 'z'
            picks=zeros(0,5);
        otherwise
    end;
    
    if doLoad
        [pa,baseName,ex]=fileparts(d.rlnMicrographName{ind});
        matName=[baseName '.mat'];
        load(matName);
        nGroups=size(p.ccImgs,3);
        mdmul=0; % display scaling not set
        displayMode=1;
    end;
    
    if displayMode==1 % Always 1 if we haven't displayed yet
        
        % Set up scaling of the filtered image
        mdf=GaussFilt(p.md,displayFilter);
        [mds,mdmul,mdadd]=imscale(mdf,256,.001);
        imaga(mds);
        
    elseif displayMode==2
        %         nb=size(
        np=size(picks,1);
        modelImage=zeros(nd,'single');
        for i=1:np
            modelImage=modelImage+ExtractImage(picks(i,3)*p.cTemplates(:,:,picks(i,4)),...
                picks(i,1:2),nd,1);
        end;
        if b=='r'
            modelImage=mdf-modelImage;
        end;
        imaga(mdmul*modelImage+mdadd);
    end;
    
    
    if boxesOn && b ~='z'
        
        % set up to search
        stackScales=ones(nGroups,1);
        stackScales(end)=nonScale;
        
        
        ccTotImgs=p.ccImgs(:,:,1);
        nd=size(ccTotImgs);
        ccTotInds=p.ccInds(:,:,1);
        ccGroups=ones(nd,'uint8');
        for i=2:nGroups
            ccXImgs=stackScales(i)*p.ccImgs(:,:,i);
            next=ccXImgs>ccTotImgs;
            ccTotImgs(next)=ccXImgs(next);
            nextInds=p.ccInds(:,:,i);
            ccTotInds(next)=nextInds(next);
            ccGroups(next)=i;
        end;
        ccGoods=any(ccGroups==p.goods,3);
        
        msz=40; % marker size
        
        ccMax=inf;
        
        ccImgWk=ccTotImgs;
        
        picks=zeros(0,4);
        
        th1=5*th0;
        
        
        hold on;
        while ccMax>th0
            [ccMax,i,j]=max2d(ccImgWk);
            good=ccGoods(i,j);
            ccImgWk=Mask(ccImgWk,[i j],blankMasks{good+1});
            if ccMax<th1
                if good % a true particle
                    marker='gs';
                    plot(i,j,marker,'markersize',msz);
                    picks(end+1,:)=[i j double(ccMax) double(ccTotInds(i,j))];
                elseif boxesOn>1
                    marker='ro';
                    plot(i,j,marker,'markersize',4);
                end;
            end;
            
        end;
        hold off;
        
    end;
    
    
    
    % %     q=struct;
    % %     p.baseName=baseName;
    % %     p.pixAds=pixAds;
    % %     p.nd=nd;
    % %     p.md=md;
    % %     p.cnd=cnd;
    % %     p.ctfParsDs=ctfParsDs;
    % %     p.cTemplates=cTemplates;
    % %     p.borders=borders;
    % %     p.goods=goods;
    % %
    % %     p.ccImgs=reshape(ccVecs,[nd nGroups]);
    % %     p.ccInds=reshape(ccInds,[nd nGroups]);
    
    txt=num2str([th0 nonScale size(picks,1)],4);
    text(1,1,txt,'verticalalignment','bottom','horizontalalignment','left',...
        'fontsize',16,'color',[1 1 0]);
    title(matName,'interpreter','none');
    %     axis off;
    
    [ix,iy,b]=MyBusyread;
    while b==0 || b='R'
        pause(0.05);
        [ix,iy,b]=MyBusyread;
    end;
    
%     [ix,iy,b]=Myginput(1);

    
    np=size(picks,1);
    
    if any(b==['nNqzR'])

        % Construct the output star file
        pk=struct;
        pk.rlnCoordinateX=picks(:,1)*ds-p.offsets(1)-1;
        pk.rlnCoordinateY=picks(:,2)*ds-p.offsets(2)-1;
        pk.rlnAutopickFigureOfMerit=picks(:,3);
        pk.rlnClassNumber=zeros(np,1);
        pk.rlnAnglePsi=zeros(np,1);
        outStarName=[outPath p.baseName '_autopick.star'];
        WriteStarFileStruct(pk,'',outStarName);
        drawnow;
        disp(['Written: ' outStarName]);
        
        
        % Compare picks
        [nms,pk0]=ReadStarFile([p.baseName '_autopick.star'],1);
        pk0=pk0{1};
        
        imags(mdf);
        hold on;
        % plot(pk.rlnCoordinateX*ds+offsets(1)+1,rlnCoordinateY*ds+offsets(2)+1),'gs');
        plot( (pk.rlnCoordinateX+offsets(1))/ds+1,(pk.rlnCoordinateY+offsets(2))/ds+1,'gs');
        plot((pk0.rlnCoordinateX+offsets(1))/ds+1,(pk0.rlnCoordinateY+offsets(2))/ds+1,'rs');
        hold off;
        
        
        
    end;
end; % while
disp('done');


return
%%

% Code to make figure for class

pk=picks(51,:);

templ=Crop(p.cTemplates(:,:,pk(4)),nd);
cc=real(ifftn(fftn(p.md).*conj(fftn(ifftshift(templ)))));
imaga(cc*300+100); hold on; plot(663,425,'gs','markersize',40); hold off;


%%

% Demonstrate signal in noise

N=randn(nd);
Nf=GaussFilt(randn(nd),.2);
si=Crop(p.cTemplates(:,:,17),nd);
X=Nf+5*si;

imags(X);

cc=real(ifftn(fftn(X).*conj(fftn(ifftshift(si)))));
% imags(cc);

%%
% Compare picks
[nms,pk0]=ReadStarFile([p.baseName '_autopick.star'],1);
pk0=pk0{1};

imags(mdf);
hold on;
% plot(pk.rlnCoordinateX*ds+offsets(1)+1,rlnCoordinateY*ds+offsets(2)+1),'gs');
plot( (pk.rlnCoordinateX+offsets(1))/ds+1,(pk.rlnCoordinateY+offsets(2))/ds+1,'gs');
plot((pk0.rlnCoordinateX+offsets(1))/ds+1,(pk0.rlnCoordinateY+offsets(2))/ds+1,'rs');
hold off;
