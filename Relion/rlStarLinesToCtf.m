function ctf=rlStarLinesToCtf(sNames,sDats,lineInds,bVal)
% function ctf=rlStarLinesToCtf(sNames,sDats,lineInds,bVal)
% Create an array of our ctf structs based on one or more star file lines.
% Run ReadStarFile to get sNames, sDats
% bVal is the B factor inserted (0 is default).
% We ignore the CtfScaleFactor and the CtfBFactor fields.
if nargin<4
    bVal=0;
end;

% Pick up the parameters that might be special to an optics group
iOptics=0; % pointer to Optics block in star file
iParticles=0;
opt=struct;
% Find the first occurrences of '*optics' or '*micrographs'/'*particles' in sNames
for i=1:numel(sNames)
    nm=sNames{i};
    if iOptics==0 && numel(strfind(nm,'optics'))>0
        iOptics=i;
    elseif iParticles==0 && ...
        (numel(strfind(nm,'micrographs'))>0 || numel(strfind(nm,'particles'))>0)
        iParticles=i;
    end;
end;

iParticles=max(1,iParticles); % has to be at least 1
pts=sDats{iParticles};

if iOptics>0 % There is an optics block
    opt=sDats{iOptics};
    pixA=opt.rlnImagePixelSize(1);
    kV=opt.rlnVoltage(1);
else
    pixA=pts.rlnPixelSize(1);
    kV=pts.rlnVoltage;
end;
lambda=EWavelength(kV);

nLines=numel(pts.rlnMicrographName);
if any(lineInds>nLines) || any(lineInds<1) % out of bounds, return nothing
    return
end;


nl=numel(lineInds);


%     CTF parameters: an array of structs
ct=struct;
for iLine=1:nl
ct.defocus=(pts.rlnDefocusU(iLine)+pts.rlnDefocusV(iLine))/2e4;
ct.deltadef=(pts.rlnDefocusU(iLine)-pts.rlnDefocusV(iLine))/2e4;
ct.theta=pts.rlnDefocusAngle(iLine)*pi/180;
ct.Cs=GetOptField('rlnSphericalAberration',iLine);
ct.alpha=GetOptField('rlnAmplitudeContrast',iLine);
ct.B=0;
ct.pixA=pixA;
ct.lambda=lambda;

if i==1
    ctf=ct;
else
    ctf(iLine)=ct;
end;
end;

    function val=GetOptField(fieldName,ind)
        % if mic.(fieldName) doesn't exist, pick it up from the opt structure.
        if isfield(pts,fieldName)
            val=pts.(fieldName)(ind);
        else
            iGrp=pts.rlnOpticsGroup(ind);
            val=opt.(fieldName)(iGrp);
        end
    end

end
