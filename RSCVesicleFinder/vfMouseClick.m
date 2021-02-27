function [h, doUpdateDisplay]=vfMouseClick(h)
% Called from a mouse click on the display, to change add or remove a
% vesicle.  What happens depends on whether the click is near the center of
% an existing vesicle.
% Right click: change to or add a good vesicle
% Left (ctrl) click: change to or add a bad vesicle
% Center (shift) click: delete the vesicle
% the function updates these variables
%   h.markedVesicleIndex (=0 if nothing changed)
%   h.goodVesImage, badVesImage
%   h.mi.vesicle (element markedVesicleIndex, if changed).


nrgn=15;  % size of box, in display pixels, allowed for clicking on
% existing vesicle
shiftVec=-h.M2(1:2,3)';
% shiftVec=h.dsMicrograph*h.dsImageShift+h.dsMicrographShift; % shift due to cropping
% in original pixels
p=get(h.axes1,'CurrentPoint');
n=h.displaySize;
p0=[p(1,1)-1 n(2)-p(1,2)]; % corrected current point
n=size(h.origImage);
p1=p0*h.ds0-shiftVec; % Selected point in global coords
if h.debug
    disp(['Local, image coords ' num2str([p0 p1])]);
end;
b=get(gcf,'SelectionType'); % normal, alt, extend are the 3 buttons.
h.changedVesicleIndex=0;   % default
doUpdateDisplay=false;     % default
badVes=h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2);

% Assign the nature of the desired new vesicle newVes
switch b
    case 'normal'  % replace or add a new 'good' vesicle
        newVesFlag=1;
    case 'alt'     % replace or add a new 'bad' vesicle
        newVesFlag=2;
    case 'extend'  % Delete a vesicle
        newVesFlag=0;
    otherwise
        newVesFlag=1;
end;

% find the index vind of the existing vesicle oldVes, and set the value of
% the newVesFlag.
%
% % vesToSearch=h.mi.vesicle.ok(:,1);
ys=h.mi.vesicle.y;
% % ys(~vesToSearch)=inf;  % exclude non-vesicles from search
dists=sqrt((h.mi.vesicle.x-p1(1)).^2 + (ys-p1(2)).^2);
oldVesFlag=0;  % default is, no vesicle found
vind=0;
v=0;

if numel(dists)>0
    [mnv, tempInd]=min(dists);
    disp(['Nearest vesicle ' num2str(mnv*h.mi.pixA,3) 'Å, index=' num2str(tempInd) ...
        '  radius ' num2str(h.mi.vesicle.r(tempInd,1)*h.mi.pixA) ' A']);
    if mnv<nrgn*h.ds0/2  % we found something within the nrgn box.
        vind=tempInd;
        oldExistFlag=h.mi.vesicle.ok(vind,1);
        oldVesFlag=1+ badVes(vind);  % 1: good vesicle; 2: bad vesicle.
    end;
end;

if oldVesFlag==0  % No existing vesicle, search the cc map
    if max(h.ccValsScaled(:))==0  % No cross correlation available
        h.markedVesicleIndex=0;  % no new vesicle
        disp('No CC map.  Find vesicles first');
        return
    end;
    %     find the nearby CC maximum and create a new vesicle entry at vind
    rctr=ceil((nrgn+1)/2);
    p1Local=round((p1+shiftVec)/h.ds0)+1;
    rgn=ExtractImage(h.ccVals,p1Local,nrgn);
    rgnAmps=ExtractImage(h.ccValsScaled,p1Local,nrgn);
    rgnRs=ExtractImage(h.ccRadii,p1Local,nrgn);
    
%     rgn=ExtractImage(h.ccValsScaled,round(p1/h.ds0)+1,nrgn);
    [mxv, i, j]=max2di(rgn);
    [~,ii,ji]=max2d(rgn);
    disp(['Amp, local coords  ' num2str([mxv i j])]);
    if (all([i j]>1) && all([i j]<nrgn) && mxv>0)  % we have a valid maximum
            amp=rgnAmps(ii,ji);
            radius=rgnRs(ii,ji);
            newcoords=h.ds0*([i j]-rctr)+p1;  % replace the position
%         Is this really a new vesicle?
        dists=sqrt((h.mi.vesicle.x-newcoords(1)).^2 + (ys-newcoords(2)).^2);
        if all(dists>h.ds0) % more than 1 pixel away
            vind=numel(h.mi.vesicle.x)+1;  % default: add a new vesicle
            h.mi.vesicle.x(vind,1)=newcoords(1);
            h.mi.vesicle.y(vind,1)=newcoords(2);
            h.mi.vesicle.s(vind,1)=amp;
%             p2=round([i j]+p1/h.ds0-rctr+1);
%             p2=max(1,min(p2,n));
%             radius=h.ccRadii(p2(1),p2(2));
            disp(['New vesicle radius ' num2str(radius*h.mi.pixA) ' A']);
            h.mi.vesicle.r(vind,1)=radius;
            h.mi.vesicle.ok(vind,1)=true;  % extend the array
            h.mi.vesicle.ok(vind,2)=(newVesFlag==1);
            % disp('Adding a vesicle: ');
            % disp(p2);
            % disp(h.mi.vesicle.s(vind))
            % disp(h.mi.vesicle.r(vind))
            v0=meMakeModelVesicles(h.mi,n,vind,0);  % create the new model
            v=real(ifftn(fftn(v0).*ifftshift(h.ctf))); % filter by ctf
            h.miChanged=1;
        end;
    end;
else  % old vesicle exists, make model and erase if on.
    if ~oldExistFlag
        oldVesFlag=0;
    end;
    v0=meMakeModelVesicles(h.mi,n,vind,0);
    v=real(ifftn(fftn(v0).*ifftshift(h.ctf))); % filter by ctf
    if (oldVesFlag==1)  % remove from the good image
        h.goodVesImage=h.goodVesImage-v;
        h.mi.vesicle.ok(vind,1:2)=false;  % mark it empty
    elseif (oldVesFlag==2)
        h.badVesImage=h.badVesImage-v;
        h.mi.vesicle.ok(vind,1:2)=false;  % mark it empty
    end;
end;

if vind>0  % a vesicle is found or created
    h.miChanged=1;
    if newVesFlag>0
        newVesFlag=mod(oldVesFlag+1,3);
    end;
    switch newVesFlag
        case 0
            disp('Erased.');
        case 1
            h.goodVesImage=h.goodVesImage+v;
            h.mi.vesicle.ok(vind,1:2)=true;
            disp('Marked good.');
        case 2
            h.badVesImage=h.badVesImage+v;
            h.mi.vesicle.ok(vind,1:2)=[true false];
            disp('Marked bad.');
    end;
end;
disp(' ');

h.markedVesicleIndex=vind;
doUpdateDisplay=true;
% disp('mc');
