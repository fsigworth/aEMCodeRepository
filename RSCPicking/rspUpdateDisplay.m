    function rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
% Version 2 handles boxes with corners more elegantly
% Here is how picks are classified
% flags        pick array index and color
% delFlag=1;      % color=1
% vesFlag=2;      % color=5
% badVesFlag=3;   % color=7;
% manFlag=16;     % color=2; **
% manRawFlag=17;  % color=6; **
% autoFlag=32;    % color=3; **
% bkdFlag=48;     % color=4
% 

if ~dis.imageMode || numel(imgs)==0
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

if dis.showGhosts > 1 % show vesicle rings 2, 3 = thin, thick lines
    x0=mi.vesicle.x/dis.ds;
    y0=mi.imageSize(2)/dis.ds-mi.vesicle.y/dis.ds;
    activeVes=all(mi.vesicle.ok(:,1:2),2);
    for i=1:numel(mi.vesicle.x)
        if activeVes(i)
            circleColor=[.1 .3 .9];
        elseif mi.vesicle.ok(i,1) % at least exists
%             circleColor=[.8 .3 .2];
            circleColor=[.8 .4 .4];
%             circleColor=[.3 .3 1];
        else
            circleColor=[.2 .2 .2];
        end;
%         lineWidth=0.7+0.3*(dis.mode>1);
        lineWidth=dis.showGhosts-1;
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
    if size(dis.classes,2)<max(ptrs)
        dis.classes=zeros(numel(ptrs),max(ptrs),'single');
    end;
    for ipt=2:numel(ptrs) % we loop through all the picks/pointers stacks except
        %         for the null stack. The index ipt sets the color also.
        coords=shiftdim(picks(ipt,1:ptrs(ipt),:),1);
        classes=shiftdim(dis.classes(ipt,1:ptrs(ipt),1),1);
        if dis.classParticlesMode && any(ipt==[2 3 6]) && numel(dis.goodClasses)>0
            clsOk=any(dis.classes(ipt,1:ptrs(ipt))'==dis.goodClasses,2);
            %             disp([sum(clsOk) numel(clsOk)]);
            nk=2;
        else
            clsOk=true(ptrs(ipt),1);
            nk=1;
        end;
        for k=1:nk % show good and bad classes
            if k==1
                colors=dis.boxColors(ipt,:);
                labelColor=dis.boxLabelColor;
                c0=coords(clsOk,:);
                cl0=classes(clsOk);
            else
                colors=min(1,dis.boxColors(ipt,:)*.7+[.5 0 0]);
                labelColor=min(1,dis.boxLabelColor*.8+repmat([.4 0 0],size(dis.boxLabelColor,1),1));
                c0=coords(~clsOk,:);
                cl0=classes(~clsOk);
            end;
            if ipt==3 %  Special autopicked case, distinguish rso/iso particles
                nj=2;  % number of box-corner types
                rsos=(c0(:,4)>0)&(c0(:,7)>0);
                c{1}=c0(~rsos,:); % the rest are iso particles.
                cl{1}=cl0(~rsos);
                c{2}=c0(rsos,:);
                cl{2}=cl0(rsos);
            else
                nj=1;
                c{1}=c0;
                cl{1}=cl0;
            end;
            
            for j=1:nj % Loop over corner types
                ourCoords=c{j};
                cls=cl{j};
                
                [bX,bY,tX,tY]=rspMakeBoxes(mi,dis,ourCoords,dis.corners(ipt,j));
                plot(bX,bY,'color',colors,'linewidth',dis.lineWidth,...
                    'HitTest','off');
                %                    Get the box 'radius' in pixels
                bPix=double(dis.currentBoxSize/mi.pixA*dis.size(1)/mi.imageSize(1));
                
                if dis.showBoxes>=3  % Show classes if present
                    for ib=1:size(ourCoords,1)
                        if cls(ib)>0
                            text(tX(ib)+bPix,tY(ib)+bPix/2,sprintf('\\bf %g',cls(ib)),...
                                'verticalAlignment','middle','horizontalAlignment','left','fontsize',dis.labelFontSize,'color',labelColor(3,:));
                        end;
                    end;
                end;
                if dis.showBoxes==4 % show picker parameters too
                    for ib=1:size(ourCoords,1)
                        if ourCoords(ib,5)>0 % if a mxCC value is given
                            %                    Show the particle amplitude on top
                            text(tX(ib),tY(ib),sprintf('\\bf%5.2g',ourCoords(ib,5)),...
                                'verticalAlignment','bottom','fontsize',dis.labelFontSize,'color',labelColor(1,:));
                            %                    Show the spectrum value (corrected by pars(11)=sFactor) on the bottom.
                            text(tX(ib),tY(ib)+bPix,sprintf('\\bf%5.2g',ourCoords(ib,8)/dis.pars(11)),...
                                'verticalAlignment','top','fontsize',dis.labelFontSize,'color',labelColor(2,:));
                        end;
                    end; % for ib
                end; % if dis.showBoxes==4
            end; % for j
        end; % for k
    end; % for ipt
end;

% Show statistics at the top of the display
labels=['  ' num2str(numel(mi.vesicle.x)) ' ' dis.imgLabels{dis.mode} ...
    '  ' num2str(dis.pars(1),2) '  ' num2str(dis.pars(10),2) ...
    ':    ' num2str(rspCountGoodParticles(picks,dis)) ' particles.    Defocus ' num2str(mi.ctf(1).defocus,3)];
if isfield(mi.ctf(1),'phi') && mi.ctf(1).phi~=0
    labels=[labels '   Phi ' num2str(mi.ctf(1).phi*180/pi,3)];
end;
text(1,1,labels,'verticalAlignment','top',...
     'fontsize',16,'fontweight','bold','color',[0 1 0]);
drawnow;
