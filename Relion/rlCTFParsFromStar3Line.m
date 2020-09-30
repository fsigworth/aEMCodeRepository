function [pars,blockIndices,nLines]=rlCTFParsFromStar3Line(sNames,sDats,iLine,B)
% function [pars,blockIndices]=rlCTFParsFromStar3Line(sNames,sDat,iLine,B)
% Handles relion version 3.1 optics groups too. Returns the indices of optics and
% micrograph/particle blocks in the sDat cell array.
% Assumes we've called [sNames,sDats]=ReadStarFile() to get the star file
% blocks. We look at the star file line iLine and returns a struct pars
% with these fields:
%  pixA lambda, defocus, Cs, B, alpha, deltadef, theta.
% pars is returned as [] if there is an error, e.g. iLine out of range.
% iLine=0 (default) returns empty pars, but assigns blockIndices and
% nLines.
% To pick up other micrograph parameters, from 
% sDats{blockIndices(2)}.rlnMicrographName{iLine}

if nargin<4
    B=50; % default value
end;

% Default outputs
pars=[];
blockIndices=[0 0];

% Pick up the parameters that might be special to an optics group
iOptics=0; % pointer to Optics block in star file
iMicrographs=0;

% Find the first occurrences of '*optics' or '*micrographs'/'*particles' in sNames
for i=1:numel(sNames)
    nm=sNames{i};
    if iOptics==0 && numel(strfind(nm,'optics'))>0
        iOptics=i;
    elseif iMicrographs==0 && ...
        (numel(strfind(nm,'micrographs'))>0 || numel(strfind(nm,'particles'))>0)
        iMicrographs=i;
    end;
end;
iMicrographs=max(1,iMicrographs); % has to be at least 1
blockIndices=[iOptics iMicrograph];

opt=struct;
if iOptics>0 % There is an optics block
    opt=sDats{iOptics};
end;
mic=sDats{iMicrographs};

nLines=numel(mic.rlnMicrographName);
if iLine>nLines || iLine<1 % out of bounds, return nothing
    return
end;


% pick up or calculate the pixel size
if isfield(mic,'rlnMicrographPixelSize') || isfield(opt,'rlnMicrographPixelSize') 
    pars.pixA=GetOptField('rlnMicrographPixelSize',iLine);
elseif isfield(opt,'rlnMicrographOriginalPixelSize')
    pars.pixA=GetOptField('rlnMicrographOriginalPixelSize',iLine);
else % compute from old-fashioned star file
    pars.pixA=1e4*mic.rlnDetectorPixelSize(iLine)/mic.rlnMagnification(iLine);
end;

%     CTF parameters
pars=struct;
pars.defocus=(mic.rlnDefocusU(iLine)+mic.rlnDefocusV(iLine))/2e4;
pars.deltadef=(mic.rlnDefocusU(iLine)-mic.rlnDefocusV(iLine))/2e4;
pars.theta=mic.rlnDefocusAngle(iLine)*pi/180;
pars.Cs=GetOptField('rlnSphericalAberration',iLine);
pars.alpha=GetOptField('rlnAmplitudeContrast',iLine);
pars.B=B;
pars.lambda=EWavelength(GetOptField('rlnVoltage',iLine));

return


    function val=GetOptField(fieldName,ind)
        % if mic.(fieldName) doesn't exist, pick it up from the opt structure.
        if isfield(mic,fieldName)
            val=mic.(fieldName)(ind);
        else
            iGrp=mic.rlnOpticsGroup(ind);
            val=opt.(fieldName)(iGrp);
        end
    end

end
