function [coords, ovMask, mxCC2]=rspAutoParticleFinder(mi,rscc,dis,mask)
% % We search for peaks in rscc.mxCC and check rscc.mxVars for excessive
% % variance, geometry, etc.
% % msub is rscc.m0-rscc.mVes, the subtracted raw image.
% % rscc.MxVesInds gives the vesicle number for each particle;
% % rscc.mxTempInds is used to get the template number.
% % each row of coords contains [x y flag vesIndex amp templ 0 0]
% % where flag = 32 (autopicked particle)
% rsoOffset (=pars(4) is how far (in angstroms) the particle center may lie
% inside the membrane center and still be counted as a bona fide right side
% out particle.   rsoOffset=0
% indicates that we take both rso and iso
% mask=1 where we reject particles in the field.
doCheckVesicleNumbers=1; % Check for consistency with rscc
% 
% np=5;
% patchMask=logical(fuzzymask(np,2,np/2,0));
% 
% zeroCC=1e-3;  % effectively zero value to search in CC
% % vsz=numel(rscc.vList)/nt;
% % vlst=reshape(rscc.vList,nt,vsz);
% % eigs=reshape(rscc.eigenImgs,nd^2,nt);
% defocus=mi.ctf(1).defocus;
% pars=dis.pars;
% % pars(12) is tFactor, the factor modifying the amplitude threshold
%   pars(11) is the sFactor, multiplies the spectrum threshold
% minAmp=max(1e-6,pars(1)+defSlope*(defocus-2))*pars(12);
minAmp=max(1e-6,dis.pars(1)*dis.pars(12)); % pars(12) is the amp scaleup factor.
% maxAmp=pars(2)*pars(12);
% maxAmp=pars(2);
% maxVar=pars(3);
% rsoOffset=pars(4)/mi.pixA;  % offset of particle center from tip (for RSO selection)
% partRadius=dis.pars(5)/mi.pixA; % blanking radius around particle, in A
overlapRadiusA=dis.pars(6); % "forbidden zone" for overlap around each vesicle, in A
% maxBob=pars(7)/mi.pixA; % maximum distance of particle center outside the vesicle.
% border=pars(8)/mi.pixA;   % border in A
% % Get the spectrum threshold (max value)
% % sThresh=pars(10)+spectSlope*(defocus-2);  % spectrum threshold
% 
% % pars(11) is tFactor, the defocus-dependent factor modifying the spectrum
% % threshold.
% sThresh=pars(10)*pars(11);  % threshold * threshFactor.
% 
n=size(rscc.mxCC,1);
dsc=mi.imageSize(1)/n(1);
% ds1=mi.imageSize(1)/n;
% n0=mi.imageSize(:)';
% angList=rscc.angleList;
% angSz=size(angList);
% angList=reshape(angList,prod(angSz(1:3)),angSz(4));
% 

% nb=2*ceil(partRadius/ds1)+4;  % size of the blanking mask box
% rscc.blankMask=1-fuzzymask(nb,2,partRadius/ds1,1);

nFound=0;

if max(rscc.mxVesInds(:))>numel(mi.vesicle.x)  % out of bounds vesicle numbers
    if doCheckVesicleNumbers
        warning('Autopicking aborted: mismatch between vesicle picker and preprocessing');
        coords=[];
        mxCC2=zeros(n,'single');
        ovMask=zeros(n,'single');
        return
    else
        rscc.mxVesInds(rscc.mxVesInds>numel(mi.vesicle.x))=0;
    end;
end;
% the amplitude values will be on the order of 1.

% Re-initialize the masked CC field
if dis.useRawAmplitudes
    mxCC2=rscc.mxCCU*dis.ccuScale;
else
    mxCC2=rscc.mxCC;
end;

% Draw a disc of r+overlapRadius around each vesicle center.  Blank the CC
% wherever these discs overlap.
[mxCC2, ovMask]=rscBlankOverlaps(mi,mxCC2,overlapRadiusA);
mxCC2=mxCC2.*(~mask);


% ---------Search for cc peaks----------
[amp, ix, iy]=max2d(mxCC2);
rscc.mxCC2=mxCC2; % this is how it's passed to the particle checker
coords=[];
while amp>minAmp
    inXY=([ix iy]-1)*dsc+1; % convert to micrograph coordinates
     [coord1,rscc.mxCC2,flags]=rspParticleChecker(amp,inXY,mi,rscc,dis);
%     [coord1,rscc.mxCC2,flags,vals,labels]=rspParticleChecker(amp,inXY,mi,rscc,dis);
%     fprintf(1,'%7.2g',coord1);
%     fprintf(1,'\n');
%     fprintf(1,' %6.3g',flags);
%     fprintf(1,'\n');
%     fprintf(1,' %6.3g',vals);
%     fprintf(1,'\n');
%     for i=1:numel(labels)
%     fprintf(1,' %s',labels{i});
%     end;
%     fprintf(1,'\n');
% pause
    if all(flags)
        nFound=nFound+1;
        if nFound==1
            coords=coord1;
            coords(1000,1)=0;
        else
            coords(nFound,:)=coord1;
        end;
    end;
%     mxCC2=Mask(mxCC2,[ix iy],blankMask);
    [amp, ix, iy]=max2d(rscc.mxCC2);
end;
coords=coords(1:nFound,:); % truncate the array

% endCC=mxCC2;

return
    
%     
%     flags(iMinAmp)=true;
%     patch=ExtractImage(rscc.mxCC2,[ix iy],np); % Copy the neighborhood
%     if all(patch(patchMask)>zeroCC) % no neighborhood points are blanked    iv=rscc.mxVesInds(ix,iy);
%         flags(iNoBlank)=true;
%         ccVar=rscc.mxVars(ix,iy);
%         if (amp < maxAmp) && (ccVar < maxVar)
%             flags.maxAmpVar=true;
%             [~, xi, yi]=max2di(rscc.mxCC2);  % get the interpolated values
%             partCtr=([xi yi]-1)*ds1;
%             ampU=rscc.mxCCU(ix,iy);  % un-normalized amplitude, goes into coords(9)
%             iv=rscc.mxVesInds(ix,iy);
%             if ~dis.spMode && iv>0 % There is a matching vesicle.
%                 flags.vesExists=true;
%                 vesCtr=[mi.vesicle.x(iv) mi.vesicle.y(iv)];
%                 vesR=mi.vesicle.r(iv);
%                 rPart=sqrt(sum((partCtr-vesCtr).^2)); % distance of particle to ves ctr
%                 %             Geometry check:
%                 %               -particle center is beyond r-rsoOffset, or else rsoOffset=0
%                 %               -particle center is within r+maxBob
%                 %               -particle center is beyond the border
%                 if (rsoOffset==0 || rPart > vesR-rsoOffset)...
%                         && (rPart <= vesR+maxBob)...
%                         && all(abs(partCtr-n0/2) < (n0/2-border))
%                     templ=single(rscc.mxTemplInds(ix,iy));
%                     rso=single(rscc.mxRsos(ix,iy));
%                     if iv==0 % no vesicle, we just write alpha
%                         rso=round(angList(max(1,templ),1)); % round to 1 degree
%                         if rso>180
%                             rso=rso-360;
%                         end;
%                     end;
%                     qp.xi=xi;
%                     qp.yi=yi;
%                     qp.mxCC2=rscc.mxCC2;
% %                     Check the residual spectrum value before inserting
% %                     the pick.
%                     [sval,xvals]=rspResidualSpectrum(mi,rscc,dis,partCtr,0);
%                     if sval<sThresh
%                         nFound=nFound+1;
% %                         put all the parameters into a row of the coords
%                         coords(nFound,10+numel(xvals))=0; % force the array bigger
%                         coords(nFound,:)=[ partCtr 32 single(iv) amp templ ...
%                                            rso sval ampU 0 xvals];
% %                                           1 2     3    4        5    6   
% %                                           7    8    9  10  +
%                     end;
%                 end;
%             elseif dis.spMode
%                 templ=single(rscc.mxTemplInds(ix,iy));
%                 alpha=round(angList(max(1,templ),1)); % round to 1 degree
%                 if alpha>180
%                     alpha=alpha-360;
%                 end;
%                 qp.xi=xi;
%                 qp.yi=yi;
%                 qp.mxCC2=rscc.mxCC2;
%                 [sval,xvals]=rspResidualSpectrum(mi,rscc,dis,partCtr,0);
%                 if sval<sThresh
%                     nFound=nFound+1;
%                     coords(nFound,:)=[ partCtr 32 single(iv) amp templ alpha sval ampU 0 xvals];
%                 end;
%                 
%             end;
%         end;
%     end;
%     % Blank the vicinity of the found peak
%     rscc.mxCC2=Mask(rscc.mxCC2,[ix iy],blankMask);
%     [amp, ix, iy]=max2d(rscc.mxCC2);
% end;
% % display the last residual
% if nFound>0
%     rspResidualSpectrum(mi,rscc,dis,coords(nFound,:),dis.showSpectrumInfo);
% end;
% coords=coords(1:nFound,:);
% endCC=rscc.mxCC2;
