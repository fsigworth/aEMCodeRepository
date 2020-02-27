function ctfPars=rlStarToCTFPars(s,i)
% Assigns values from the ith row of a star file (i.e. from the ith element of
% each of the struct fields s as returned by ReadStarFile) to our standard
% CTF parameters struct. 

pars.B=50;
pars.ampFactor=1;
pars.phaseShift=0;
fieldNames=fieldnames(s);

ctfPars=struct;
ctfPars.defocus= (s.rlnDefocusU(i) + s.rlnDefocusV(i))/2e4;
ctfPars.deltadef=(s.rlnDefocusU(i) - s.rlnDefocusV(i))/2e4;
ctfPars.theta=s.rlnDefocusAngle(i)*pi/180;
ctfPars.phi=sget('rlnPhaseShift',pars.phaseShift)*pi/180;
% ctfPars.phi=s.rlnPhaseShift(i)*pi/180;
ctfPars.alpha=sget('rlnAmplitudeContrast');
% ctfPars.alpha=s.rlnAmplitudeContrast(i);
ctfPars.lambda=EWavelength(s.rlnVoltage(i));
ctfPars.Cs=s.rlnSphericalAberration(i);
ctfPars.resLimit=sget('rlnCtfMaxResolution');
ctfPars.estB=sget('rlnCtfBfactor',pars.B);
ctfPars.ccc=s.rlnCtfFigureOfMerit(i);

ctfPars.B=pars.B;
ctfPars.ampFactor=1;

ctfPars.pixA=s.rlnDetectorPixelSize(i)/s.rlnMagnification(i)*1e4;

% delete unassigned fields
pNames=fieldnames(ctfPars);
for j=1:numel(pNames)
    if numel(ctfPars.(pNames{j}))==0
        ctfPars=rmfield(ctfPars,pNames{j});
    end;
end;


% local function

function x=sget(name,x) % x argument is default value.
    if nargin<2
        x=[];
    end;
    if any(strcmp(name,fieldNames))
        x=s.(name)(i);
    end;    
end
end
