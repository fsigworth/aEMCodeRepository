function miOut=rsMergeVesicleList2(mi,origMi)
% Merge the original vesicle list with the newly found ones in mi.
% This allows the particle assignments to vesicles to remain the same when
% the vesicle finder is run after particle picking.
% Copy the particle picks from origMi.  Update the vesicle list, keeping
% the same index for the vesicles that match the old vesicle positions.

maxOffset=100;  % angstroms, cutoff for identifying the "same" vesicle.
updateOldVesicles=1;  % write new parameters for particle-associated vesicles

% Test if any particles have been picked.
npv=0;
if isfield(origMi.particle,'picks') && size(origMi.particle.picks,1)>0
    npv=sum(origMi.particle.picks(:,4)>0);  % number of pointers to vesicles
end;

% Make a boolean array of replaced vesicle params
nvOld=numel(origMi.vesicle.x);
nvNew=numel(mi.vesicle.x);
if nvOld<1 || nvNew<1
    miOut=mi;
    return
end;
origMi.vesicle.ok(nvOld,4)=false; % force the size
mi.vesicle.ok(nvNew,4)=false;
freeOld=true(nvOld,1);  % Mark unused slots in the original mi file
freeNew=true(nvNew,1);
% hasParticle=false(nvOld,1);

miOut=mi;   % By default we copy the new information, but we'll replace
%           particles and some vesicles.

if npv>0  % There are pointers to vesicles, so we adjust the vesicle entries 
    disp(['Updating ' num2str(npv) ' marked vesicles.']);
    miOut.vesicle=origMi.vesicle;  % By default, copy the old vesicle info
%     Remove higher-order vesicle radius and amplitudes
    miOut.vesicle.s(:,2:end)=[];
    miOut.vesicle.r(:,2:end)=[];
    miOut.vesicle.ok=(false(size(miOut.vesicle.ok))); % mark them all bad.
    miOut.particle=origMi.particle; % preserve the particles.
    newxs=mi.vesicle.x(:);
    newys=mi.vesicle.y(:);
    maxErr=maxOffset/mi.pixA;
% size(miOut.vesicle.ok)
% size(mi.vesicle.ok)
    for oldVesIndex=1:nvOld  % scan through all the old vesicles, and
%                               search for corresp. new ones.
        oldx=origMi.vesicle.x(oldVesIndex);
        oldy=origMi.vesicle.y(oldVesIndex);
        oldr=origMi.vesicle.r(oldVesIndex);
        %             Find the closest new one.
        dists=(newxs-oldx).^2+(newys-oldy).^2;
        [minDist, newVesIndex]=min(dists);
        if minDist<maxErr^2 % we found a new one close enough
            %                 Check the radius too
            if abs(mi.vesicle.r(newVesIndex)...
                    -origMi.vesicle.r(oldVesIndex)) > maxErr
                warning(['Vesicle radii don''t match: '...
                    num2str([mi.vesicle.r(newVesIndex) oldr])]);
            end;
            if updateOldVesicles
                % copy the relevant new vesicle's info into the old position
                % of miOut.
                miOut.vesicle.x(oldVesIndex)=mi.vesicle.x(newVesIndex);
                miOut.vesicle.y(oldVesIndex)=mi.vesicle.y(newVesIndex);
                miOut.vesicle.r(oldVesIndex,1)=mi.vesicle.r(newVesIndex,1);
                miOut.vesicle.s(oldVesIndex,1)=mi.vesicle.s(newVesIndex,1);
                miOut.vesicle.ok(oldVesIndex,:)=mi.vesicle.ok(newVesIndex,:);
            end;
            %                 Mark this vesicle index as used
            freeOld(oldVesIndex)=0;  % This position is used
            freeNew(newVesIndex)=0;  % Won't need to copy this one.
            newxs(newVesIndex)=inf;  % Don't allow this one to be found again
        end;  % if minDist
    end; % for oldVesIndex

    numNewToInsert=sum(freeNew);  % See how many new vesicles still to be copied.
    numOldFree=sum(freeOld);      % How many unassigned positions
    if numOldFree<numNewToInsert
        %         Extend the list of vesicles so there will be enough for
        %         copies of the unassigned new ones.
        freeOld=[freeOld; true(numNewToInsert-numOldFree,1)];
    end;
    oldIndices=find(freeOld);
    newIndices=find(freeNew);
    oldIndices=oldIndices(1:numel(newIndices));
        
    %     Copy all the new params into the unused slots.
    miOut.vesicle.x(oldIndices)=mi.vesicle.x(newIndices);
    miOut.vesicle.y(oldIndices)=mi.vesicle.y(newIndices);
    miOut.vesicle.r(oldIndices,1)=mi.vesicle.r(newIndices,1);
    miOut.vesicle.s(oldIndices,1)=mi.vesicle.s(newIndices,1);
    miOut.vesicle.ok(oldIndices,:)=mi.vesicle.ok(newIndices,:);
end;  % if there are particles
