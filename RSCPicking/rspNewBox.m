function [picks, ptrs, rscc, mi, doUpdate]=rspNewBox(mi,rscc,dis,picks,ptrs,coords,clickType) % insert a manual pick
% Interpret a click. coords are padded micrograph coordinates, 0-based.
% The doUpdate flag is set to 1 when a bad vesicle is marked, so any
% deleted particles are also marked.

% flags        pick array index and color
delFlag=1;      % color=1
vesFlag=2;      % color=5
badVesFlag=3;   % color=7;
manFlag=16;     % color=2
manRawFlag=17;  % color=6;
% autoFlag=32;  % color=3 (not used here)
bkdFlag=48;     % color=4
minAmp=1e-3;    % amp value that is essentially zero.

requiresPreprocessorText='SimpleRSP:ChangedFlags';

flags=0;
flagLabels=[];
vals=zeros(1,9,'single');
nsets=numel(ptrs);

dsc=dis.ds;  % Downsample factor of mxCC

doUpdate=0;

coords(3)=0;
if numel(coords)<10
    coords(10)=0;   % The written coords are: [x y type vesInd cc template rso alpha]
end;

allInds=zeros(nsets,1);
minDistances=zeros(nsets,1);
minDistances(1)=inf;
allInds(1)=1;
% We find the minimum distances (in original pixels) of the current coords
% from each of the picks stacks.
for i=2:nsets  % we ignore i=1 which are deleted particles
    dists=hypot(picks(i,1:ptrs(i),1)-coords(1),...
        picks(i,1:ptrs(i),2)-coords(2));
    if numel(dists)<1
        dists=inf;
    else
        dists(picks(i,:,3)==0)=inf;  % Don't count bad vesicle coords
    end;
    [minDistances(i), allInds(i)]=min(dists);
end;

[globalMin, setInd]=min(minDistances);

if isfield(mi,'vesicle') && isfield(mi.vesicle,'x') && numel(mi.vesicle.x)>0
    vesDists=hypot(mi.vesicle.x-coords(1), mi.vesicle.y-coords(2));
    [~, vesInd]=min(vesDists);
else
    vesInd=0;
end;




switch clickType
    case {1,'.'}  % simple click, should be a manual pick, goes into picks(2,..)
        %         We'll select it if we're not at an old location
        if globalMin<dis.minDist/mi.pixA; % minDist is in angstroms, convert to pixels.
            DoErase;
        else
            %                  if it's beyond an existing pick, assume it's a new one.
            coords(3)=manRawFlag;  % mark as manually picked
            if clickType==1
                %             Check against the CC map
                nrgn=ceil(0.3*dis.pars(20)/(dsc*mi.pixA))*2+1;  % compute from box size:
                %                                       an odd number, 2/3 of box size
                rctr=ceil((nrgn+1)/2);
                %                 disp(['Starting coords ' num2str(coords(1:2))]);
                localCoords=round(coords(1:2)/dsc)+1;
                rgn=ExtractImage(rscc.mxCC2,localCoords,nrgn);
                %                 Take the maximum over the region.
                [mxv, i, j]=max2di(rgn);
                if ~(any([i j]<1) || any([i j]>nrgn) || mxv<=minAmp)  % we have a valid maximum
                    inXY=dsc*([i j]-rctr)+dsc*(localCoords-1);  % refined micrograph position
                    %                 disp(['  CC max coords ' num2str(inXY)]);
                    
                    %                     xL=i-rctr+localCoords(1);
                    %                     yL=j-rctr+localCoords(2);
                    %                     Figure out whether this would be an autopick
                    %                     origCoords=coords;
                    %          [xL yL]
                    [coords,mxCC2,flags,vals,flagLabels]=rspParticleChecker(mxv,inXY,mi,rscc,dis);
                    %                 disp(['   Final coords ' num2str(coords(1:2))]);
                    coords
                    %                 disp(['corrected coords: ' num2str([origCoors(1:2) ' -> ' newcoords])]);
                    if flags(5)>0  % it's associated with a vesicle
                        coords(3)=manFlag;
                    else
                        coords(3)=manRawFlag;
                    end;
                end;
                if coords(3)==manFlag  % This has corrected coordinates
                    boxColor=2;
                    ptrs(2)=ptrs(2)+1;
                    picks(2,ptrs(2),1:numel(coords))=coords;
                else
                    coords(4)=vesInd;
                    boxColor=6;
                    ptrs(6)=ptrs(6)+1;
                    picks(6,ptrs(6),1:numel(coords))=coords;
                end;
                ShowInfo('New: ',coords,flags,vals,flagLabels);
                iiso=((coords(7)>0) & (coords(4)>0))+1; % always zero if no vesicle.
                [bX,bY,tX,tY]=rspMakeBoxes(mi,dis,coords,dis.corners(boxColor,iiso));
                plot(bX,bY,'color',dis.boxColors(boxColor,:),'linewidth',dis.lineWidth);  % draw a box.
                if dis.showBoxes>2
                    text(tX,tY,sprintf('\\bf%5.2g',coords(5)),'verticalAlignment','bottom','fontsize',12,'color','g');
                    bPix=double(dis.currentBoxSize/dis.pixA);
                    text(tX,tY+bPix,sprintf('\\bf%5.2g',coords(8)),'verticalAlignment','top','fontsize',12,'color','y');
                end;
            end;
        end;
    case {2,'g'}  % shift-click or center button or 'k' erases.  This goes into picks(1,...)
        DoErase;
        
    case 3  % ctrl-click marks a non-particle. picks(4,...)
        coords(3)=bkdFlag;  % mark as a background particle
        ptrs(4)=ptrs(4)+1;
        picks(4,ptrs(4),1:numel(coords))=coords;
        [bX, bY]=rspMakeBoxes(mi,dis,coords,dis.corners(4));
        plot(bX,bY,'color',dis.boxColors(4,:));  % draw a box.
        
    case 'W'  % get liposome information
        vdist=inf;
        if numel(mi.vesicle.x)>0 % vesicles exist
            vdist=hypot(mi.vesicle.x(:)-coords(1),mi.vesicle.y(:)-coords(2));
        end;
        [minVes, indVes]=min(vdist);
        if minVes<dis.minDist    % We are close to a vesicle center
            %             List the vesicle info
            disp(['Vesicle index: ' num2str(indVes)]);
            disp(['Radius (Å): ' num2str(mi.vesicle.r(indVes,1)*mi.pixA)]);
            disp('Terms:');
            i=3;
            while i<size(mi.vesicle.r,2) && mi.vesicle.r(indVes,i)~=0
                disp([num2str(i-1) '  ' num2str(abs(mi.vesicle.r(indVes,i))*mi.pixA)]);
                i=i+1;
            end;
            %             ptrs(5)=ptrs(5)+1;
            %             boxColor=5;                    %             Insert the vesicle index
            %             coords(1)=mi.vesicle.x(indVes);
            %             coords(2)=mi.vesicle.y(indVes);
            %             coords(3)=vesFlag;
            %             coords(4)=indVes;
            %             picks(5,ptrs(5),:)=coords;
            %             [bX, bY]=rspMakeBoxes(mi,dis,coords,dis.corners(5));
            %             plot(bX,bY,'color',dis.boxColors(boxColor,:));  % draw a box.
        end;
    case {'v' 'x'}  % toggle good/bad vesicle
        doUpdate=1;  % We'll let the calling program update the display.
        vdist=inf;
        if numel(mi.vesicle.x)>0 % vesicles exist
            vdist=hypot(mi.vesicle.x(:)-coords(1),mi.vesicle.y(:)-coords(2));
        end;
        [minVes, indVes]=min(vdist);
        %         [minVes indVes]
        if minVes<2*dis.minDist    % We are close to a vesicle center
            
            if clickType=='v' % was a bad vesicle, mark it good.
                mi.vesicle.ok(indVes,2)=~mi.vesicle.ok(indVes,2);
                disp(mi.vesicle.ok(indVes,:));
                %Make a log entry to say that we need an update
                logIndex=numel(mi.log); % First check if the last entry is already this message
                if ~strncmp(mi.log{logIndex},requiresPreprocessorText,numel(requiresPreprocessorText))
                    logIndex=numel(mi.log)+1;
                end;
                %                 Either way, update the log entry with the time stamp.
                mi.log{logIndex,1}=[requiresPreprocessorText ' ' TimeStamp];
                
            else % 'x' typed.
                mi.vesicle.ok(indVes,2)=0; % mark the vesicle bad
                ptrs(7)=ptrs(7)+1;   % make a vesicle entry
                coords(1)=mi.vesicle.x(indVes);
                coords(2)=mi.vesicle.y(indVes);
                coords(3)=badVesFlag;
                coords(4)=indVes;
                picks(7,ptrs(7),1:numel(coords))=coords;
                [picks, ptrs]=rspDeleteBadVesiclePicks(picks,ptrs);
                %                 if dis.mode==8 % we're looking at vesicle models
                %                     mVesDel=meMakeModelVesicles(mi,dis.ndis,indVes);
                %                     disp(mi.vesicle.ok(indVes,:));
                %                     vesPresent=mi.vesicle.ok(indVes,1);
                %                     if vesPresent
                %                         rscc.mVes=rscc.mVes-mVesDel;
                %                     else
                %                         rscc.mVes=rscc.mVes+mVesDel;
                %                         DoErase;
                %                     end;
                %                     mi.vesicle.ok(indVes,1)=~vesPresent;
                %                 end;
            end;
        end;
        
        %         else % draw a miniature box to show that we've responded.
        %             origBoxSize=dis.pars(20);
        %             dis.pars(20)=dis.pars(20)/2;
        %             [bX, bY]=rspMakeBoxes(mi,dis,coords);
        %             dis.pars(20)=origBoxSize;
        %             plot(bX,bY,'color',dis.boxColors(7,:));  % draw a box.
        %          end;
        
    case 'w'  % get info about a particle        if globalMin<dis.minDist  % We clicked on an existing marker
        if globalMin<dis.minDist  % We clicked on an existing marker
            partInd=allInds(setInd);  % Look up the particle info
            coords=squeeze(picks(setInd,partInd,:))';  % Mark it deleted.
            ShowInfo('Info: ',coords);
        end;
end;  % switch

    function DoErase
        if globalMin<dis.minDist*2  % We clicked near an existing marker; allow 2x larger error.
            %             allInds(setInd)
            partInd=allInds(setInd);
            picks(setInd,partInd,3)=0;  % Mark it deleted.
            boxLoc=picks(setInd,partInd,:);
            ShowInfo('Deleted: ',boxLoc);
            boxLoc(3)=1;  % mark it erased
            [bX, bY]=rspMakeBoxes(mi,dis,boxLoc);
            plot(bX,bY,'color',dis.boxColors(1,:));  % draw a black box.
            if setInd==3  % special case:  it was an auto-picked particle.
                %                     We need to store the location.
                ptrs(1)=ptrs(1)+1;  % we'll make an entry in the deleted stack
                picks(1,ptrs(1),:)=picks(setInd,partInd,:);
                picks(1,ptrs(1),3)=delFlag;  % Mark the entry as deleted.
            end;
        end;
    end

    function ShowInfo(str,coords,flags,vals,flagLabels)
        if any(coords(1:6)<0)
            return
        end;
        if nargin<5
            vals=0;
        end;
        coords=coords(:);
        coords(10)=0; % extend the array
        
        if dis.listParticleInfo  % if we've turned this on
            rcoords=round(coords);
            fprintf(1,'\n');
            fprintf(1,[str 'x= %d  y= %d   flag= %d ves= %d\n      cc=%5.2f ref= %d  rso= %d  spec= %d\n'],...
                rcoords(1),rcoords(2),rcoords(3),rcoords(4),coords(5),rcoords(6),rcoords(7),rcoords(8));
            if any(vals ~= 0) % we tested the particle
                fmtStr=repmat(' %6d',1,numel(flags));
                fprintf(1,[fmtStr '\n'],flags);
                fmtStr=repmat(' %6.4g',1,numel(vals));
                fprintf(1,[fmtStr '\n'],vals);
                for i=1:numel(flagLabels)
                    fprintf(1,'%s ',flagLabels{i});
                end;
                fprintf(1,'\n');
            end;
        end;
        if dis.showSpectrumInfo
            rspResidualSpectrum(mi,rscc,dis,coords,1);
        end;
        
        
    end % showInfo
end % main function
