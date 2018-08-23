% rlCheckDefBest.m
% run after rlCheckDefocus to find repreesentative particles
% Used to create images for BPS 2017 poster
% Used to create jpegs of picked particles in images
writeJpegs=0;
drawBoxes=1;

jpegDir='/Users/fred/Documents/ ppt stuff/My talks/2017/';

%%%%%%
cd('/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/Stack2');
load sq81_mining.mat
disp('Loading the si file');
load sq81_1p192m2+tsi.mat
%%%%%%

cd('/Users/fred/EMWork/Hideki/160909/KvLipo121_2retrack');
mergeDir='Merged/';
microDir='Micrograph_p/';
iOffset=1223;

% cd('/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/');
% mergeDir='MergedZ/';
% iOffset=0;

subText={'mz.tif' 'mvz.tif'};
flagText={'m' 's' 'e1' 'e2'};

% micros=1:10000;
% for i=micros
%     j=i+iOffset;
%     mi=si.mi{j};
%     disp([num2str([i mnParts(j) mCls(j,2)],4) '   ' mi.baseFilename]);
% end;

%%
figure(5);

        dis.org=[0 0];
        dis.currentBoxSize=180;
        dis.size=[1920 1920];
        bColorErased=[.3 .3 .3];
        bColorMan=[1 .7 0];  % yellow
        bColorManNotRefined=[1 .55 .3];
        bColorAuto=[.8 1 0];  % greenish
%         bColorAut2=[.6 .8 .4];
        bColorBkd=[.4 1 .8]; % cyan
        bColorVesicle=[.2 0 1]; % blue vesicle
        bColorBadVesicle=[1 .1 .3]; % red bad vesicle

        dis.boxColors=[bColorManNotRefined; bColorAuto; bColorMan;
            bColorMan; bColorVesicle; bColorManNotRefined; bColorBadVesicle];
        dis.boxLabelColor=[0 1 0; .5 1 0];
        dis.corners=[1  1  1 1  .41 1.2 .41
                     0 .8 .7 1  .41 1.2 .41]'; %  'eroded' box corners
        dis.lineWidth=2;
%         dis.lineWidth=0;
        

% index=input('Index? ');
for index=33:33  % micrograph index
    j=index+iOffset;
    mi=si.mi{j};
    cParts=cParticles{j};
    sParts=pInds(cParts);  % si particle indices
    mParts=si.miParticle(sParts);  % mi particle indices
    tnParts=size(mi.particle.picks,1);

    ncls=max(cls);
    clsParts=false(tnParts,ncls);
    for i=1:ncls
        clsParts(mParts,i)=cls(cParts)==i;
    end;    
    
    baseName=[mergeDir mi.baseFilename];
    %     micName=mi.imageFilenames{1};
    %     m1=ReadEMFile([microDir micName]);
    for k=1
        if k<3
        nm=[baseName subText{k}];
        else
            nm=[microDir mi.imageFilenames{k-2}];
        end;
        m=ReadEMFile(nm);
        %     mv=ReadEMFile([baseName 'mvz.tif']);
        
        % Make the figure
%         figure(5);
        clf;
        
        set(gcf,'toolbar','none');
        % ,'resize','off');
        
        % Main display
        pos=get(gcf,'Position');
        xsiz=pos(3);
        ysiz=pos(4);
        
        xsiz=960;
        ysiz=960;
        
        set(gcf,'Position',[pos(1:2) xsiz ysiz]);
        axsiz=xsiz;
        aysiz=ysiz;
        dis.ax1=axes('units','pixels','position',[2 3 axsiz aysiz],'ticklength',[0 0]);
        
        dis.ds=mi.imageSize(1)/xsiz;
        
        
        % Model vesicle code
        
        % v=meMakeModelVesicles(mi,dis.size,0,0,0);
        
        
        % mf=GaussFilt(Downsample(m-mv,dis.size),.1*dis.ds);
        % mul=754;
        % add=128;
        % ms=(mf*mul+add);
        
        mf=GaussFilt(Downsample(m,xsiz),.08*dis.ds);
        if numel(mf)<2
            break
        end;
        % [ms,mul,add]=imscale(mf,256,.0001);
        if any(k==[1 3 4])
%             [ms,mul,add]=ImscaleQuadrant(mf,256,.0001,3);
            [ms,mul,add]=imscale(mf,256,.0001);
        else  % re-use the scaling.
            ms=mf*mul+add;
        end;
        
        imaga(ms);
        
        % imwrite(uint8(rot90(ms)),[jpegDir 'M1Img.jpg']);
        % Add boxes to the figure
        axis off
        title(mi.baseFilename,'interpreter','none');
        %
coords=mi.particle.picks;
        if drawBoxes && size(coords,1)>0
            
            hold on;
            labelNos=zeros(size(coords,1));
            % goodParts=cls(cParts);
            for colorIndex=1:ncls
                    drawnParts=clsParts(:,colorIndex);
                for rso=0:1
                    coordInds=find(any(clsParts,2));
                    labelNos(coordInds)=sParts;
                    coords(:,3)=coords(:,7)==rso & drawnParts;  % rso only
                    [bX,bY,tX,tY]=rspMakeBoxes(mi,dis,coords,dis.corners(3,rso+1));
                    plot(bX,dis.size(2)-bY,'color',dis.boxColors(colorIndex,:),'linewidth',dis.lineWidth,...
                        'HitTest','off');
                    for ip=find(drawnParts)
%                         disp([double(coords(ip,1)/dis.ds) double(ysiz-coords(ip,2)/dis.ds)]);
%                         plot(double(coords(ip,1)/dis.ds),double(ysiz-coords(ip,2)/dis.ds),'ws');
                      text(double(coords(ip,1)/dis.ds),double(coords(ip,2)/dis.ds-dis.currentBoxSize/(2*dis.ds*mi.pixA)),...
                          num2str(labelNos(ip)),'color',dis.boxColors(colorIndex,:),...
                          'fontsize',12,'horizontalalignment','center','verticalalignment','top');
                    end;
                end;
            end;
            hold off;
        end;
if writeJpegs
    outDir='/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/Stack2/';
        figDir='Figures/';
        sn=sprintf('%04d%s_',j,flagText{k});
        outName=[sn mi.baseFilename '.jpg'];
        disp(outName);
        print([outDir figDir outName],'-djpeg','-r0');
end;
    end;  % for k
end;