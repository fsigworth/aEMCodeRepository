% NISMapAnalysis.m

doLoad=0;

addpath ~/aEMCodeRepository/aLibs/HealpixLib
cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia')
if doLoad
    disp('Loading...');
    load NISMapData5.mat
    disp('done.');
end;
% loads: m5v p5r ptrsI ligandLabels dsv nv s names
%%
figSizes = ...
    [      0           0           0           0
    -2517         795        1120         736
    -1685          66        1104         734
    -860         897         560         420
    -2905          78        1538         755
    -2625         891         855         650 ];
%    % set figure sizes
%    for i=2:size(figSizes,1)
%        figure(i);
%        set(gcf,'Position',figSizes(i,:));
%    end;
%    figure(1);
%
mode='radial';
% mode='ligands';
showAllAngles=1;
doSave=1;
doOptimizeShifts=0; % choose the cdIShifts to maximize density at ion positions.

% ptrsL=ptrsL1;
yLimsAll=[-1.5 2.5];
yLims1=[-.5 1.6];

%        mv=circshift(m4v,[0 0 9]); % need to shift to match peaks.
%         [0 0 9] for J322 works a bit better.
mv=m5v;
        mapText=names.pdb5;
        [~,modelName]=fileparts(mapText);
        figPrefix=modelName;
        if numel(figPrefix)>17
            figPrefix=figPrefix(1:17);
        end;
        if doOptimizeShifts
            figPrefix=[figPrefix '_sh'];
        end;
        pr=p5r;
        ptrsI=ptrsI5;
        nIons=numel(ptrsI);
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
    NISMapOptimizeShifts;
    disp('Ion position shifts')
    disp(cdIShifts/dsv);
end;
% return

cdIs=cdIs+cdIShifts;


% cdLs=[pr.X(ptrsL)' pr.Y(ptrsL)' pr.Z(ptrsL)']*dsv+ctrv;

healPixOrder=2;

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
        
        if showAllAngles
            hpEulers = HealpixGenerateSampling(healPixOrder, 'scoord');
            vecs = SphToCart(hpEulers);
            nDirs=size(vecs,1)/2;
            vecs=vecs(1:nDirs,:); % pick up non-reflected points.
            lines=zeros(2*range+1,nDirs);
        end;
        
        for j=1:nIons
            ptr=ptrsI(j);
            if pr.flipped(ptr)
                flipMark='*';
            else
%                 flipMark=''; % bizarre error that causes problem with title.
                flipMark=' ';
            end;
            ionLabel=[pr.element{ptr} ' chain ' pr.chainID{ptr} ' ' ...
                num2str(pr.resNum(ptr)) flipMark];
            ionLabel(15)=' ';
            ionLabels{j}=ionLabel;
            mysubplot(2,nIons,j);
            cdI=cdIs(j,:);
            
            imags(mv(:,:,round(cdI(3)))); % show the section
            hold on;
            plot(cdI(1),cdI(2),'yo','markersize',20);
            hold off;
            axis off;
            ionLabelX=[ionLabel '  [' num2str(cdIShifts(j,:)/dsv,2) ']'];
            ionLabelX
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
            [rMean,rMedian,rVals,rVals]=Radial3(mLoc,[]);
            rVals=rMedian;
            yVals(1:range,j)=rVals(end:-1:2);
            yVals(range+1:2*range+1,j)=rVals;
            
            if j==3
                lw=1.5;
            else
                lw=1;
            end;
            
            if showAllAngles
                plot(xVals,yVals(:,j),'k-','linewidth',lw);
                hold on;
                for k=1:nDirs
                    lines(:,k)=ExtractLine3(mv,nPts,vecs(k,:),cdI);
                end;
                plot(xVals,lines);
                yLims=yLimsAll;
            else
                plot(xVals,yVals(:,j),'-','color',colors(j,:),'linewidth',lw);
                hold on;
                yLims=yLims1;
            end;
            plot(xVals,0*xVals,'k-');
            hold off;
            grid on;
            axis([-inf inf yLims]);
            ylabel('Map density');
            xlabel('Radius, Å');
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
        
        figure(4);
        set(gcf,'position',figSizes(4,:));
        clf;
        title('Means');
        plot(xVals,yVals,'linewidth',1);
        hold on;
        plot(xVals,0*xVals,'k-');
        hold off;
        axis([min(xVals) max(xVals) yLims]);
        xlabel('Radial distance from center, Å');
        ylabel('Map density');
        grid on;
        legTxt=ionLabels;
        legend(legTxt);
        if doSave
            print('-djpeg','-r300',[figPrefix '_Radial.jpg']);
                        disp(['Wrote ' figPrefix '_Radial.jpg']);

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
