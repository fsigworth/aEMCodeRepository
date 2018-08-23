function [picks, ptrs]=rspDeleteBadAutoPicks(dis,picks,ptrs)
% scan each deleted box and delete from auto picks.
% picks(1,...) are deleted boxes
% picks(3,...) are auto picks

if ptrs(3)<1  % no auto-picked particles to scan
    return
end;

for i=1:ptrs(1)  % scan all deleted coordinates.
    coords=picks(1,i,:);
%     coords(:)'
    autoDist=hypot(picks(3,1:ptrs(3),1)-coords(1),...
        picks(3,1:ptrs(3),2)-coords(2));
    [minAuto indAuto]=min(autoDist);
    if minAuto<dis.minDist
        % Draw a black box in the old location
%         [bX bY]=rspMakeBoxes(mi,dis,picks(3,indAuto,:),dis.boxRadius);
%         plot(bX,bY,'color',dis.boxColors(1,:));  % draw a black box.
          %     Mark the auto-pick as deleted.
        picks(3,indAuto,3)=0;
    end;
end;
