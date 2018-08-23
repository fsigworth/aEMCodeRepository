function rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
% Version 2 handles boxes with corners more elegantly
if ~dis.imageMode
    return
end;

% Draw the initial image
img=imgs(:,:,dis.mode);
% img=dis.imageScaleUp*(imgs(:,:,dis.mode)-120);

cimg=rspMakeRGB(dis,img,masks);
hold off;
hi=image(cimg);
set(hi,'HitTest','off');
% GetClick('init');
hold on;

if dis.showBoxes==4 % show vesicles too
    x0=mi.vesicle.x/dis.ds;
    y0=mi.imageSize(2)/dis.ds-mi.vesicle.y/dis.ds;
    activeVes=all(mi.vesicle.ok,2);
    for i=1:numel(mi.vesicle.x)
        if activeVes(i)
            circleColor=[.1 .3 .9];
        else
%             circleColor=[.8 .3 .2];
            circleColor=[.8 .4 .4];
        end;
        lineWidth=0.7+0.3*(dis.mode>1);
            [xs,ys]=CircleLineSegments(mi.vesicle.r(i,:)/dis.ds,4);
            plot(xs+x0(i),ys+y0(i),'-','color',circleColor,'lineWidth',lineWidth);
%             maxAmp=max(VesicleCircAmps(mi.vesicle.s(i,:,1),64));
%             str=sprintf('%5.2g [%4.2g %4.2g]',i,maxAmp*1e4,mi.vesicle.s(i,1,1)*1e4);
%             text(double(x0(i)),double(y0(i)+min(ys)),str,'color',[1 1 0],...
%                 'verticalAlignment','bottom','horizontalAlignment','center','fontsize',12);
    end;
end;

% coords( x y flag vesInd mxCC templInd rso )

if dis.showBoxes
    for ipt=2:numel(ptrs) % we loop through all the picks/pointers stacks except
        %         for the null stack. The index i0 sets the color also.
        coords=shiftdim(picks(ipt,1:ptrs(ipt),:),1);
        %             default coordinates and number of box variations
        c{1}=coords;
        nj=1; % number of box variations
        if ipt==3 %  Special autopicked case, distinguish rso/iso particles
            nj=2;  % number of box-corner types
            %           we'll separate coords into 2 parts
            %           Rso particles: if the vesicle is assigned, we take only the coords(7) flagged
            %           particles.
            rsos=(coords(:,4)>0)&(coords(:,7)>0);
            c{1}=coords(~rsos,:); % the rest are iso particles.
            c{2}=coords(rsos,:);
        end;

        for j=1:nj % Draw all the boxes of category ipt
            [bX,bY,tX,tY]=rspMakeBoxes(mi,dis,c{j},dis.corners(ipt,j));
            plot(bX,bY,'color',dis.boxColors(ipt,:),'linewidth',dis.lineWidth,...
                'HitTest','off');
            if dis.showBoxes==3  % Show text with boxes
                %                    Get the box 'radius' in pixels
                bPix=double(dis.currentBoxSize/mi.pixA*dis.size(1)/mi.imageSize(1));
                for ib=1:size(coords,1)
                    if coords(ib,5)>0 % if a mxCC value is given
                        %                    Show the particle amplitude on top
                        text(tX(ib),tY(ib),sprintf('\\bf%5.2g',coords(ib,5)),...
                            'verticalAlignment','bottom','fontsize',12,'color',dis.boxLabelColor(1,:));
                        %                    Show the spectrum value (corrected by pars(11)=sFactor) on the bottom.
                        text(tX(ib),tY(ib)+bPix,sprintf('\\bf%5.2g',coords(ib,8)/dis.pars(11)),...
                            'verticalAlignment','top','fontsize',12,'color',dis.boxLabelColor(2,:));
                    end;
                end;
            end;
        end;
        
    end;
end;

% Show statistics at the top of the display
labels=['  ' dis.imgLabels{dis.mode} ...
    '  ' num2str(dis.pars(1),2) '  ' num2str(dis.pars(10),2) ...
    ':    ' num2str(ptrs(3)) ' particles.    Defocus ' num2str(mi.ctf(1).defocus,3)];
text(1,1,labels,'verticalAlignment','top',...
     'fontsize',16,'fontweight','bold','color',[1 .8 .9]);
drawnow;
