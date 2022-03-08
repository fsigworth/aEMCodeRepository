% NISMapAnalysisGen.m
%  Version for making figures 17-feb-22
% Call this after running NISMapLoaderGen.m

doLoad=1;

addpath ~/aEMCodeRepository/aLibs/HealpixLib
cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia/')
dataDir='220216/';

% ----- Here we pick which dataset we're working on, I- or ReO- -----
% dataName='NISMapDataI.mat';
% isReo=0;
% maxNIons=inf;
% yLimsAll=[-1.2 2.8];
% yLimsMean=[-.7 1.8];

dataName='NISMapDataREO.mat';
isReo=1;
maxNIons=3;
maxNIonsInMeanPlot=maxNIons;
yLimsAll=[-1.5 5];
yLimsMean=[-1 5];
% -----------------------------------------------------
fileName=[dataDir dataName];
disp(fileName);

if doLoad
    disp(['Loading ' fileName '...']);
    load(fileName);
    disp('done.');
end;
% loads: m5v p5r ptrsI ligandLabels dsv nv s names sites



%%
figSizes = ...
    [      0           0           0           0
    -1517         795        1120         736
    -1355          66        1104         734
    -600         897         560         420
    -1905          78        1538         755
    -1625         891         855         650 ];
%    % set figure sizes
%    for i=2:size(figSizes,1)
%        figure(i);
%        set(gcf,'Position',figSizes(i,:));
%    end;
%    figure(1);
%
mode='radial';
% mode='ligands';
medianMode=0; % Use the median rather than mean over directions.
showAllAngles=1;
showSites=1;

doSave=1;

doOptimizeShifts=1; % choose the cdIShifts to maximize density at ion positions.
iVersion=4;

mv=m5v;
mapText=names.pdb;
[~,modelName]=fileparts(mapText);
figPrefix=modelName;
if numel(figPrefix)>17
    figPrefix=figPrefix(1:17);
end;
if doOptimizeShifts
    figPrefix=[figPrefix '_sh'];
end;
if ~medianMode
    figPrefix=[figPrefix '_me'];
end;
figPrefix=[dataDir figPrefix num2str(iVersion)];
pr=p5r;
nIons=min(maxNIons,numel(ptrsI));
ctrv=ceil((nv+1)/2);

% % % simulated map for checking cc of map to model
% fakeM=zeros(nv,nv,nv,'single');
% for i=1:numel(pr.X)
%     fakeM(round(pr.X(i)*dsv+ctrv),round(pr.Y(i)*dsv+ctrV),round(pr.Z(i)*dsv+ctrV))=1;
% end;
% fakeM=GaussFilt(fakeM,.1);
% cc=fftshift(real(ifftn(fftn(circshift(mv,[0 0 0])).*conj(fftn(fakeM)))));
% ShowSections(cc);
% [val,ishift]=max3di(cc)
%


% convert from zero-based original voxels to pixel indices in mv

cdIs=[pr.X(ptrsI)' pr.Y(ptrsI)' pr.Z(ptrsI)']*dsv+ctrv; % padded pixels

% cdIs(2,:)=cdIs(2,:)+[1 0 0];
% cdIs(3,:)=cdIs(3,:)+[0 0 -1]; % shift Y by 1/4 voxel.

cdIShifts=zeros(nIons,3);

if doOptimizeShifts
    disp('Optimizing shifts.');
    NISMapOptimizeShifts;
    if isReo
        cdIShifts(3,:)=0; % no shifts for the fake Na2
    end;
    disp('Ion position shifts, A')
    for i=1:nIons
        disp([num2str(cdIShifts(i,:)/dsv) '   ... ' ligandLabels5{i}]);
    end;
end;
% return

cdIs=cdIs+cdIShifts;

% Get the C-alpha norms
cIndices=find(strncmp('C',ligandLabels5,1));
nCs=numel(cIndices);
pkVals=zeros(nCs,1);
for i=1:nCs
    cdI=round(cdIs(cIndices(i),:));
    pkVals(i)=mv(cdI(1),cdI(2),cdI(3));
end;
disp(['C-alpha peak values before norm: ' num2str(pkVals',3)]);
pkNorm=mean(pkVals)/2.2; % Makes Na+ to have density around 1.
mv=mv/pkNorm;
% cdIs(cIndices,:)=[]; % suppress these peaks.
nIonsDis=min(nIons,size(cdIs,1)); % now show only these ions.


% cdLs=[pr.X(ptrsL)' pr.Y(ptrsL)' pr.Z(ptrsL)']*dsv+ctrv;

healPixOrder=2;
%%
switch mode
    case 'radial'
        figure(3-showAllAngles);
        set(gcf,'position',figSizes(3-showAllAngles,:));
        clf;
        range=5*dsv; % range, in pixels
        nx=2*range+3; % makes the radial function of size range+1
        xVals=(-range:range)*s.pixA/dsv; % actual angstroms
        nPts=2*range+1;
        yVals=zeros(nPts,nIons);
        ionLabels=cell(nIons,1);


        for j=1:nIons
            ptr=ptrsI(j);
            if pr.flipped(ptr)
                flipMark='*';
            else
                %                 flipMark=''; % bizarre error that causes problem with title.
                flipMark=' ';
            end;

            if showAllAngles
                hpEulers = HealpixGenerateSampling(healPixOrder, 'scoord');
                vecs = SphToCart(hpEulers);
                nDirs=size(vecs,1)/2;
                vecs=vecs(1:nDirs,:); % pick up non-reflected points.
                sel=1:2:nDirs;
                vecs=vecs(sel,:); % half as many
                nDirs=nDirs/2;
                lines=zeros(2*range+1,nDirs);
            end;

            if showSites % Show the coordinating atom sites
                cdSs=[pr.X(sites(j).ptrs)' pr.Y(sites(j).ptrs)' pr.Z(sites(j).ptrs)']*dsv+ctrv;
                nDirs=size(cdSs,1);
                vecs=cdSs-repmat(cdIs(j,:),nDirs,1); % vectors from sites to ion
                lines=zeros(2*range+1,nDirs);
            end;
            ionLabel=[pr.element{ptr} ' chain ' pr.chainID{ptr} ' ' ...
                num2str(pr.resNum(ptr)) flipMark];
            ionLabel(15)=' ';
            ionLabel(ionLabel==char(0))=' ';
            ionLabels{j}=ionLabel;
            mysubplot(2,nIons,j);
            %             subplot(2,nIons,j);
            cdI=cdIs(j,:);

            imags(mv(:,:,round(cdI(3)))); % show the section
            hold on;
            plot(cdI(1),cdI(2),'yo','markersize',20);
            hold off;
            axis off;
            shiftStr=num2str(cdIShifts(j,:)/dsv,2);
            ionLabelX=[ionLabel '  ' shiftStr ];
            %             titleText=[figPrefix '   ' ionLabelX]
            %             ionLabelX=[ionLabel '  [' num2str(cdIShifts(1,:)/dsv,2) ']'];
            %             ionLabelX
            if j==1
                title([figPrefix '   ' ionLabelX],'interpreter','none');
                %                 title(ionLabelX,'interpreter','none');
            else
                title(ionLabelX,'interpreter','none');
            end;

            mysubplot(2,nIons,nIons+j)
            colors=get(gca,'colororder');

            %             sVecs=rand(10,3)-.5;
            %             for i=1:size(sVecs,1)
            mLoc=ExtractVolumeInterp(mv,cdI,nx);
            [rMean,rMedian]=Radial3(mLoc,[]);
            if medianMode
                rVals=rMedian;
            else
                rVals=rMean;
            end;
            yVals(1:range,j)=rVals(end:-1:2);
            yVals(range+1:2*range+1,j)=rVals;


            if showAllAngles
                textXs=zeros(nDirs,1);
                for k=1:nDirs
                    %                     if sites.dirs(k)
                    lines(:,k)=ExtractLine3(mv,nPts,vecs(k,:),cdI);
                    textX=1;
                    %                     else
                    %                         lines(:,k)=ExtractLine3(mv,nPts,-vecs(k,:),cdI);
                    %                         textX=-1;
                    %                     end;
                    %                         if showSites
                    %                             %                         lines(:,sites.dirs==0)=-lines(:,sites.dirs==0);
                    %                             ctr=floor(nPts/2+1);
                    %                             dx=round(nPts*.35);
                    %                             text(textX*4.5,lines(ctr+dx*textX,k),sites.res{k},...
                    %                                 'HorizontalAlignment','center','verticalAlignment','middle');
                    %                         end;
                end;
                lw1=1.2;
                plot(xVals,lines,'-','linewidth',lw1); % and superimpose the colored individual lines
                hold on;
                yLims=yLimsAll;
            else
                plot(xVals,yVals(:,j),'-','color',colors(j,:),'linewidth',lw1);
                hold on;
                yLims=yLims1;
            end;
            %        Plot the mean
                lwMean=3;
                plot(xVals,yVals(:,j),'-','color',[.6 .6 .6],'linewidth',lwMean); % plot the mean/median
            plot(xVals,0*xVals,'k-'); % baseline
            hold off;
            if showSites
                legend([sites(j).res {'mean' ''}],'location','northwest');
            end;
%             plot(xVals,0*xVals,'k-'); % baseline
%             hold off; 

            grid on;
            axis([-inf inf yLims]);
            ylabel('Map density');
%             xlabel('Radius, Å');
            xlabel('Distance toward ligand, Å');
        end; % for j
        if doSave
            if showAllAngles
                figTxt='_Radial_multi.jpg';
            else
                figTxt='_Radial1.jpg';
            end;
            print('-djpeg','-r300',[figPrefix figTxt]);
            disp(['Wrote ' figPrefix figTxt]);
        end;
% ----------------- Show the means in one plot--------------
        figure(4);
        lw2=2;
        nIonsDis=min(nIons,maxNIonsInMeanPlot);
        set(gcf,'position',figSizes(4,:));
        clf;
        title('Means');
        plot(xVals,yVals(:,1:nIonsDis),'linewidth',lw2);
        hold on;
        plot(xVals,0*xVals,'k-');
        hold off;
        axis([min(xVals) max(xVals) yLimsMean]);
        xlabel('Radial distance from center, Å');
        ylabel('Map density');
        grid on;
        legTxt=ionLabels;
        legend(legTxt);
        if doSave
            jpegName=[figPrefix '_Radial.jpg'];
            print('-djpeg','-r300',jpegName);
            disp(['Wrote ' jpegName]);

        end;


    case 'ligands'
        %         show where we are
        % figure(10);
        % imags(m1vcr(:,:,cds(3)));
        % hold on;
        % plot(pr.X(ptrsL)*us*vs+ctrv,pr.Y(ptrsL)*us*vs+ctrv,'go','markersize',10);
        % plot(cds0(1),cds0(2),'yo','markersize',10);
        % hold off;

        nL=numel(ptrsL);
        range=7.5*dsv; % range, in pixels

        nPts=2*range+1;
        xVals=(-range:range)*s.pixA/dsv;
        lines=zeros(nPts,nL);
        points=zeros(nL,1);
        figure(5);
        set(gcf,'position',figSizes(5,:));
        clf;
        colors=get(gca,'colororder');
        for k=1:nL
            jIon=ligandIons(k);
            cdI=cdIs(jIon,:);
            cdL=cdLs(k,:);

            % % %          Testing the ligand location.
            %             cdLr=round(cdL);
            %             mv(cdLr(1),cdLr(2),cdLr(3))=2;
            % % %

            vec=cdL-cdI;
            %             line increases in the direction of the ligand
            lines(:,k)=ExtractLine3(mv,nPts,vec,cdI);
            points(k)=range+1+round(sqrt(vec*vec'));
            %              plot(xVals(points(k)),lines(points(k),k),'k+','markersize',10);

            mysubplot(2,nL,k);
            imags(mv(:,:,round(cdL(3))));
            axis off;
            hold on;
            plot(cdL(1),cdL(2),'yo','markersize',20);
            plot(cdI(1),cdI(2),'w+','markersize',10);
            hold off;
            if k>1
                title(ligandLabels(k));
            else
                title([mapText '   ' ligandLabels{1}],'interpreter','none');
            end;
            mysubplot(2,nL,k+nL);
            plot(xVals,lines(:,k),'-','linewidth',1,'color',colors(k,:));
            hold on;
            plot(xVals(points(k)),lines(points(k),k),'.','markersize',20,'color',colors(k,:));
            plot(xVals,0*xVals,'k-');
            hold off;
            grid on;
        end;


        figure(6);
        set(gcf,'position',figSizes(6,:));
        clf;
        ax=gca;
        plot(xVals,lines,'-','linewidth',1);
        %         plot(xVals,lines);
        %
        hold on;
        for k=1:nL
            plot(xVals(points(k)),lines(points(k),k),'.','markersize',20,'color',colors(k,:));
        end;
        plot(ax.XLim,[0 0],'k-');
        plot([0 0],ax.YLim,'k-');
        hold off;
        grid on;
        ligLabels=ligandLabels;
        ligLabels(end+1:end+nL+1)={' '};
        legend(ligLabels);
        % We assume wer'e dealisng iwth only one ion.!!
        ionLabel=[pr.element{ptrsI(jIon)} ' chain ' pr.chainID{ptrsI(jIon)}];
        title([mapText '  ' ionLabel],'interpreter','none');

    otherwise
        disp(['Unrecognized mode: ' mode]);
end;

return;

%% Pick up figure sizes
sizes=zeros(6,4);
for i=2:6
    figure(i);
    sizes(i,:)=get(gcf,'Position');
end;

%
%% Show the vicinity of the NaF4 spot
saveIt=1;
j=2;
cdI=cdIs(j,:);
[mvAbs,mulr,addr]=imscale(mv(:,:,round(cdI(3))),256,1e-4);
nax=19 % number of z-sections to show
na0=(nax+1)/2; % offset of center
nextr=7; % no. of points in extracted subimage
xtoffs=(nextr+1)/2-3;
for iOs=1:nax
    mysubplot(4,5,iOs);
    mvd=mv(:,:,round(cdI(3))+iOs-na0);
    mvs=ExtractImage(mvd,round(cdI(1:2)+3),nextr);
    [maxv,ix,iy]=max2di(mvs);
    imaga(mvd*mulr+addr); % show the section
    hold on;
    plot(cdI(1),cdI(2),'yo','markersize',20);
    plot(cdI(1)+ix-xtoffs,cdI(2)+iy-xtoffs,'ro','markersize',15);
    hold off;
    axis off;
    title(['z= ' num2str(iOs-na0) '  peak val= ' num2str(maxv,3)]);
end;
if saveIt
    print('-djpeg','-r300',[figPrefix '_NaF4_slices.jpg']);
end;
