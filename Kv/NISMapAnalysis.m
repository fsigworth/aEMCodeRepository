% NISMapAnalysis.m

mode='radial';


cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia')
load NISMapData
% loads: m1v m2v p1r ptrsI ptrsL ligandLabels dsv nv s names
%%

mv=m1v; % Choose the Na_I map
nIons=numel(ptrsI);

figure(2);

ctrv=ceil((nv+1)/2);

switch mode
    case 'radial'
        range=5*dsv;
        nx=4*range+3;
        xVals=(-range:range)*s.pixA/dsv;
        yLims=[0 1.5];
%         Convert ion coordinates to max pixels
        cdIs=[p1r.X(ptrsI)' p1r.Y(ptrsI)' p1r.Z(ptrsI)']*dsv/s.pixA+ctrv; 
        for j=1:nIons
            ptr=ptrsI(j);
            ionLabel=[p1.element{ptr} ' chain ' p1.chainID{ptr}];
            mysubplot(2,nIons,j);
            cds=cdIs(j,:);
            imags(mv(:,:,round(cds(3)))); % show the section
            hold on;
            plot(cds(1),cds(2),'yo','markersize',20);
            hold off;
            axis off;
            title(ionLabel);
            %             nx=NextNiceNumber(2*range+5)
            
% %             Tweak the ion location
%             if j>1 % the first one is problematic
%                 dLoc=ExtractVolumeInterp(m1vcr,cds,7);
%                 [val, dcds]=max3di(dLoc);
%                 cds1=cds+dcds-4; % tweaked cds0
%                 disp([j cds1-cds0]);
%             else
%                 cds1=cds0;
%             end;
            mLoc=ExtractVolumeInterp(mv,cds,nx);
%             Tune up the center coordinate
%             imags(sum(mLoc,3));
            [rMean,rMedian,yVals,rVals]=Radial3(mLoc,[]);
            
%             % Correct the r=0 point
%             rMean(1)=0*rMean(1)+1*mean(yVals{2}(1:6)); % 1/2 times r=1 points
%             rMean(2)=4/3*rMean(2)-1/3*mean(yVals{2}(1:6)); % remove the r=1 points
%             
%             sMean=[rMean(range+1:-1:1); rMean(2:range+1)];
%             sMedian=[rMedian(range+1:-1:1); rMedian(2:range+1)];
%             lines(:,1,1,j)=sMean;
%             lines(:,2,1,j)=sMedian;
%             legTxt{j}=ionLabel;
%             %             radLine=zeros(2*range+1);
%             %             radLine(range+1:end)=rMed(1:range);
%             %             radLine(range:-1:1)=rMed(2:range-1);
%             %
            mysubplot(2,nIons,nIons+j)
            plot(xVals,[rMean rMedian]);
            plot(xVals,rMedian,'linewidth',1);
            hold on
            plot(xVals,0*xVals,'k-');
            hold off;
            grid on;
            axis([-inf inf yLims]);
            ylabel('Map density');
            xlabel('Radius, Å');
        end; % for j
%     case 'ligands'
%         if j~=ligandIon  % only want one ion.
%             continue;
%         end;
% %         show where we are
% % figure(10);
% % imags(m1vcr(:,:,cds(3)));
% % hold on;
% % plot(p1r.X(ptrsL)*us*vs+ctrv,p1r.Y(ptrsL)*us*vs+ctrv,'go','markersize',10);
% % plot(cds0(1),cds0(2),'yo','markersize',10);
% % hold off;
% 
%         figure(bigFigure+4);
% %         plot(xVals,0*xVals,'k-');
% %         hold on;
%         nPts=2*range+1;
%         xVals=(-range:range)*s.pixA/(us*vs);
%         lines=zeros(nPts,nL);
%         points=zeros(nL,1);
% %         colors=get(gca,'colororder');
%         cds0
%         for k=1:nL
%             ptrL=ptrsL(k);
%             cdl=[p1r.X(ptrL) p1r.Y(ptrL) p1r.Z(ptrL)]*us*vs+ctrv
%             vec=cdl-cds0;
%             lines(:,k)=ExtractLine3(m1vcr,nPts,vec,cds0);
%             points(k)=range+1-round(sqrt(vec*vec'));
% %              plot(xVals(points(k)),lines(points(k),k),'k+','markersize',10);
%         end;
%             plot(xVals,lines,'-','linewidth',1);
% %         plot(xVals,lines);
% %         
%         hold on;
%         for k=1:nL
%                      plot(xVals(points(k)),lines(points(k),k),'k+','markersize',10);
%         end;
%         hold off;
%             grid on;
%         legend([ligandLabels' {'' '' ''}]);
%             
%     case 'angs'            
%             rngs0=[{0} {0} {0}];
%             for k=1:3 % x,y,z trace
%                 rngsh=rngs0;
%                 rngsh(k)={-range:range};
%                 lines(:,k,i,j)=squeeze(m1vcr(cds(1)+rngsh{1},cds(2)+rngsh{2},cds(3)+rngsh{3}));
%             end;
%             %         mysubplot(2,nIons,j);
%             %         linesToPlot=reshape(lines,2*range+1,3*nAngs,nIons);
%             %         plot( (-range:range)*s.pixA/(us*vs) , linesToPlot(:,ltp,j) ,'linewidth',1);
%             %         axis([-inf inf min(lines(:)) max(lines(:))]);
%             %         legend(legtxt);
%             % Show the vicinity of the ion
%             iAng=mod(i-1,2);
%             kAng=floor((i-1)/2);
%             col=2*(j-1)+iAng+1;
%             rowA=1+kAng;
%             rowB=3+kAng;
%             %         disp([i j rowB col]);
%             nr=4;
%             nc=nIons*2;
%             mysubplot(nr,nc,nc*(rowA-1)+col);
%             curves=1:3;
%             if i>1
%                 lines(:,5-i,i,j)=0; % zero out the redundant trace
%             end;
%             plot( xVals, lines(:,:,i,j) ,'linewidth',1);
%             hold on;
%             plot( xVals , 0*xVals, 'k-','linewidth',.5);
%             hold off;
%             axis([-inf inf min(lines(:)) max(lines(:))]);
%             if i==1 && j==1
%                 legend('x','y','z');
%             end;
%             grid on;
%             if i==1 % first angle
%                 title(ionLabel);
%             end;
%             
%             mysubplot(nr,nc,nc*(rowB-1)+col);
%             imags(m1vcr(:,:,cds(3)));
%             hold on;
%             plot(cds(1),cds(2),'yo','markersize',20);
%             hold off;
%             axis off;
%         end; % showRadialAverage
%     drawnow;
%     end; % for j=1:nIons
end; % for i-1:nAngs
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
