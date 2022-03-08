% NISMapOptimizeShifts.m
%  Introduce small shifts into the ion coordinates to maximize the
%  radially-averaged density.
%  Called from early in NISMapAnalysis5.m

% Variables expected:
%     s, dsv,nIons, ptrsI, cdIs, mv
%     we'll return cdIShifts(nIons,3)
doDisplay=1;

cdIShifts=zeros(nIons,3);
range=5*dsv; % range, in pixels
nx=2*range+3; % makes the radial function of size range+1
xVals=(-range:range)*s.pixA/dsv; % actual angstroms
nPts=2*range+1;
yVals=zeros(nPts,nIons);
ionLabels=cell(nIons,1);
steps=0.2*[1 1 1];
% figure(10);
for j=1:numel(ptrsI)
    ptr=ptrsI(j);
    cdI=cdIs(j,:);
    P=zeros(1,3);
    P=Simplex('init',P,steps);
    for i=1:100; % iterations
        cdI=cdIs(j,:)+P;
        mLoc=ExtractVolumeInterp(mv,cdI,nx);
        [rMean,rMedian]=Radial3(mLoc,[]);
        rVals=rMedian;
        yVals(1:range,j)=rVals(end:-1:2);
        yVals(range+1:2*range+1,j)=rVals;
        P=Simplex(-rVals(1));
        if mod(i,25)==0 && doDisplay
            disp(P);
%             plot(xVals,yVals);
%             drawnow;
        end;
    end;
    cdIShifts(j,:)=P;
%     pause(0.5);
disp('------');
end;

