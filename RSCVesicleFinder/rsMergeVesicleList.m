function mi=rsMergeVesicleList(mi,origMi)
% Reassign particles to vesicle indices.  Scan the origMi.particle.picks
% for vesicles identified with picks.  Find the closest vesicle centers in
% the new mi corresponding to those in the origMi, and change the
% origMi.particle.picks(:,4) entries which are the vesicle numbers. Then
% copy the origMi.particle.picks to mi and return it.


maxOffset=100;  % angstroms
% Test if any particles have been picked.
np=0;
if isfield(origMi.particle,'picks') && size(origMi.particle.picks,1)>0
    np=size(origMi.particle.picks,1);
    npv=sum(origMi.particle.picks(:,4)>0);  % number of identified vesicles
end;
nv1=numel(mi.vesicle.x);  % number of new vesicles
numUpdated=0;
numCopied=0;

if np>0  % There are particles
    disp(['Updating ' num2str(np) ' particle assignments to vesicles.']);
    newxs=mi.vesicle.x(:);
    newys=mi.vesicle.y(:);
    maxErr=maxOffset/mi.pixA;
    
    for i=1:np
        oldVesIndex=origMi.particle.picks(i,4);
        if oldVesIndex>0 && oldVesIndex<=numel(origMi.vesicle.x)
            oldx=origMi.vesicle.x(oldVesIndex);
            oldy=origMi.vesicle.y(oldVesIndex);
            oldr=origMi.vesicle.r(oldVesIndex);
            dists=(newxs-oldx).^2+(newys-oldy).^2;
            [minDist newVesIndex]=min(dists);
            if minDist<maxErr^2 % we found a new one close enough
%                 Check the radius too
                if abs(mi.vesicle.r(newVesIndex)...
                        -origMi.vesicle.r(oldVesIndex)) > maxErr
                    warning(['Vesicle radii don''t match: ' num2str(...
                        [mi.vesicle.r(newVesIndex)...
                        origMi.vesicle.r(oldVesIndex)])]);
                end;
%                 disp([oldVesIndex newVesIndex]);
                oldMi.particle.picks(i,4)=newVesIndex;
                numUpdated=numUpdated+1;
            else % copy the old vesicle position and radius
                numCopied=numCopied+1;
                nv1=nv1+1;
                mi.vesicle.x(nv1)=oldx;
                mi.vesicle.y(nv1)=oldy;
                mi.vesicle.r(nv1)=oldr;
                mi.vesicle.s(nv1)=median(mi.vesicle.s);
            end;
        else  % old vesicle index wasn't valid
            origMi.particle.picks(i,4)=0;  % no vesicle assigned
        end;
    end;
    numUpdated
    numCopied
    mi.particle.picks=origMi.particle.picks;  % copy the particles
end;  % if there are particles
    