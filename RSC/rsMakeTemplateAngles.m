function templateAngles=rsMakeTemplateAngles(runInfo)
% return the nangs x 3 array, with gamma varying most quickly.
angleStep=runInfo.angleStep;
nBeta=round(180/angleStep+1);
dGamma=360/(runInfo.nGamma*runInfo.symmetry);
nangs=nBeta*runInfo.nGamma;
templateAngles=zeros(nangs,3);
i=1;
for iBeta=0:nBeta-1
    for iGamma=0:runInfo.nGamma-1
       templateAngles(i,:)=[0  iBeta*angleStep iGamma*dGamma];
       i=i+1;
    end;
end;
