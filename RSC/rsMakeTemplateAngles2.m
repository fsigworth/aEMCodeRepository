function [templateAngles runInfo]=rsMakeTemplateAngles2(runInfo,gammaMode)
if nargin<2
    gammaMode='uniform';
end;

% return the nangs x 3 array, with gamma varying most quickly.
angleStep=runInfo.angleStep;
nBeta=round(180/angleStep+1);
dGamma=360/(runInfo.nGamma*runInfo.symmetry);
nangs=nBeta*runInfo.nGamma;
templateAngles=zeros(nangs,3);
i=1;
switch gammaMode
    case 'uniform'
for iBeta=0:nBeta-1
    runInfo.gammaPointers(iBeta+1,:)=[i runInfo.nGamma];
    for iGamma=0:runInfo.nGamma-1
       templateAngles(i,:)=[0  iBeta*angleStep iGamma*dGamma];
       i=i+1;
    end;
end;
        
    case 'quantized'        
        
        for iBeta=0:nBeta-1
    beta=iBeta*angleStep;
    gammaStride=min(floor(1/sind(beta)),runInfo.nGamma);
    runInfo.gammaPointers(iBeta+1,:)=[i runInfo.nGamma];

    nGamma=);
    for iGamma=0:nGamma-1
       templateAngles(i,:)=[0  iBeta*angleStep iGamma*dGamma];
       i=i+1;
    end;
end;
