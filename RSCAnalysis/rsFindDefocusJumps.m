% rsFindDefocusJumps.m

allMisName='Picking_9/allMis_holemasked_i2+f.mat';
% load(allMisName);

defs=miExtractFieldValues(allMis,'ctf.defocus');
res=miExtractFieldValues(allMis,'ctf.resLimit');
nmi=numel(defs);

cDefs=defs;
cDefs(res>7 | defs>4)=NaN;
plot(cDefs);

xDefs=cDefs;
for i=1:nmi
    if isnan(xDefs(i))
        xDefs(i)=xDefs(i-1);
    end;
end;
fDefs=GaussFilt(xDefs,.02);
dfDefs=diff(fDefs);
modelPts=round(216:484.8:nmi);
transPts=modelPts;
np=numel(modelPts);
offs=30;
for i=1:np
    search=modelPts(i)-offs:modelPts(i)+offs;
    [~,dp]=max(dfDefs(search));
    if dp==1 || dp==2*offs+1
        dp=offs+1;
    end;
    transPts(i)=modelPts(i)-offs-1+dp;
end;
model0=zeros(nim-1,1);
model0(modelPts)=.1;
model=zeros(nim-1,1);
model(transPts)=.1;
xDefs1=xDefs(2:nim);
plot([model dfDefs model0 xDefs1/10]);
% These match at the points after the defocus jump. So,
% transPts are the points just before the defocus jump.
defJumpStarts=transPts;
save Picking_9/defJumpStarts.mat defJumpStarts