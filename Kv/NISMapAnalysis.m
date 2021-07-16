% NISMapAnalysis.m

addpath ~/aEMCodeRepository/aLibs/HealpixLib
cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia')
disp('Loading...');
load NISMapData
disp('done.');
% loads: m1v m2v p1r ptrsI ptrsL ligandLabels dsv nv s names

figSizes = ...
    [      0           0           0           0
       -2517         795        1120         736
       -1685          66        1104         734
        -860         897         560         420
       -2905          78        1538         755
       -2625         891         855         650 ];
%%
mode='radial';
% mode='ligands';
showAllAngles=1;

% mv=m1v; % Choose the Na_I map
% mapText='Na_I_Map';
mv=m2v;
mapText='Na_Map';

nIons=numel(ptrsI);


ctrv=ceil((nv+1)/2);
cdIs=[p1r.X(ptrsI)' p1r.Y(ptrsI)' p1r.Z(ptrsI)']*dsv+ctrv; % padded pixels
cdIs(1,:)=cdIs(1,:)+[1 -1 0];
cdIs(2,:)=cdIs(2,:)+[0 -1 0];

cdLs=[p1r.X(ptrsL)' p1r.Y(ptrsL)' p1r.Z(ptrsL)']*dsv+ctrv;
hpOrder=2;

switch mode
    case 'radial'
        figure(3-showAllAngles);
        set(gcf,'position',figSizes(3-showAllAngles,:));
        clf;
        range=5*dsv; % range, in pixels
        nx=2*range+3; % makes the radial function of size range+1
        xVals=(-range:range)*s.pixA/dsv; % actual angstroms
        yLims=[-.5 1.6];
        nPts=2*range+1;
        yVals=zeros(nPts,nIons);
        ionLabels=cell(nIons,1);
        
        if showAllAngles
            hpEulers = HealpixGenerateSampling(hpOrder, 'scoord');
            vecs = SphToCart(hpEulers);
            nDirs=size(vecs,1)/2;
            vecs=vecs(1:nDirs,:); % pick up non-reflected points.
            lines=zeros(2*range+1,nDirs);
        end;
        
        for j=1:nIons
            ptr=ptrsI(j);
            ionLabel=[p1r.element{ptr} ' chain ' p1r.chainID{ptr}];
            ionLabels{j}=ionLabel;
            mysubplot(2,nIons,j);
            cdI=cdIs(j,:);

            imags(mv(:,:,round(cdI(3)))); % show the section
            hold on;
            plot(cdI(1),cdI(2),'yo','markersize',20);
            hold off;
            axis off;
            if j==1
                title([mapText '   ' ionLabel],'interpreter','none');
            else
                title(ionLabel,'interpreter','none');
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
                yLims=[-1.5 2.5];
            else
                plot(xVals,yVals(:,j),'-','color',colors(j,:),'linewidth',lw);
                hold on;
            end;
            plot(xVals,0*xVals,'k-');
            hold off;
            grid on;
            axis([-inf inf yLims]);
            ylabel('Map density');
            xlabel('Radius, Å');
        end; % for j

figure(4);
set(gcf,'position',figSizes(4,:));
clf;
yLims=[-.5 1.6];
title('Means');
plot(xVals,yVals,'linewidth',1);
hold on;
plot(xVals,0*xVals,'k-');
hold off;
axis([min(xVals) max(xVals) yLims]);
legTxt=ionLabels;
legend(legTxt);

    case 'ligands'
        %         show where we are
        % figure(10);
        % imags(m1vcr(:,:,cds(3)));
        % hold on;
        % plot(p1r.X(ptrsL)*us*vs+ctrv,p1r.Y(ptrsL)*us*vs+ctrv,'go','markersize',10);
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
        ionLabel=[p1r.element{ptrsI(jIon)} ' chain ' p1r.chainID{ptrsI(jIon)}];
        title([mapText '  ' ionLabel],'interpreter','none');

%     case 'angs'
%         
%         
%         for j=1:nIons
%             lines(:,k)=ExtractLine3(m1vcr,nPts,vec,cdI);
%             points(k)=range+1-round(sqrt(vec*vec'));
%             %              plot(xVals(points(k)),lines(points(k),k),'k+','markersize',10);
%             
%             mysubplot(2,nI,k);
%             imags(mv(:,:,round(cdL(3))));
%             hold on;
%             plot(cdL(1),cdL(2),'yo','markersize',20);
%             plot(cdI(1),cdI(2),'w+','markersize',10);
%             hold off;
%             title(ligandLabels(k));
%             
%             mysubplot(2,nL,k+nL);
%             plot(xVals,lines(:,k),'-','linewidth',1,'color',colors(k,:));
%             hold on;
%             plot(xVals(points(k)),lines(points(k),k),'.','markersize',20,'color',colors(k,:));
%             plot(xVals,0*xVals,'k-');
%             hold off;
%             grid on;
%         end;
%         
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



%                     rngs0=[{0} {0} {0}];
%                     for k=1:3 % x,y,z trace
%                         rngsh=rngs0;
%                         rngsh(k)={-range:range};
%                         lines(:,k,i,j)=squeeze(m1vcr(cds(1)+rngsh{1},cds(2)+rngsh{2},cds(3)+rngsh{3}));
%                     end;
%                     %         mysubplot(2,nIons,j);
%                     %         linesToPlot=reshape(lines,2*range+1,3*nAngs,nIons);
%                     %         plot( (-range:range)*s.pixA/(us*vs) , linesToPlot(:,ltp,j) ,'linewidth',1);
%                     %         axis([-inf inf min(lines(:)) max(lines(:))]);
%                     %         legend(legtxt);
%                     % Show the vicinity of the ion
%                     iAng=mod(i-1,2);
%                     kAng=floor((i-1)/2);
%                     col=2*(j-1)+iAng+1;
%                     rowA=1+kAng;
%                     rowB=3+kAng;
%                     %         disp([i j rowB col]);
%                     nr=4;
%                     nc=nIons*2;
%                     mysubplot(nr,nc,nc*(rowA-1)+col);
%                     curves=1:3;
%                     if i>1
%                         lines(:,5-i,i,j)=0; % zero out the redundant trace
%                     end;
%                     plot( xVals, lines(:,:,i,j) ,'linewidth',1);
%                     hold on;
%                     plot( xVals , 0*xVals, 'k-','linewidth',.5);
%                     hold off;
%                     axis([-inf inf min(lines(:)) max(lines(:))]);
%                     if i==1 && j==1
%                         legend('x','y','z');
%                     end;
%                     grid on;
%                     if i==1 % first angle
%                         title(ionLabel);
%                     end;
%
%                     mysubplot(nr,nc,nc*(rowB-1)+col);
%                     imags(m1vcr(:,:,cds(3)));
%                     hold on;
%                     plot(cds(1),cds(2),'yo','markersize',20);
%                     hold off;
%                     axis off;
%                 end; % showRadialAverage
%             drawnow;
%             end; % for j=1:nIons
%
% if strcmp(mode,'radial')
%     figure(2*bigFigure);
%     plot(xVals,squeeze(lines(:,1,1,:)));
%     legend(legTxt);
%     title('Means');
%     if bigFigure==2
%         figText='Na_I_Map';
%         maxY=1.5;
%     else
%         figText='Na_Map';
%         maxY=1.5;
%     end;
%     figure(3);
%     plot(xVals,squeeze(lines(:,2,1,:)),'linewidth',1);
%     hold on;
%     plot(xVals,0*xVals,'k-');
%     hold off;
%     axis([min(xVals) max(xVals) -.5 maxY]);
%
%
%     legend(legTxt);
%     title([figText ' Medians'],'interpreter','none');
%     xlabel('Radius from atom center, Å');
%     ylabel('Map density');
% end;
%
