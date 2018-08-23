

%% TestAutoParticleFinder
%
% We search for peaks in mxCC and check mxVars for excessive variance.
% MxVesInds gives the vesicle number for each particle; mxTempInds can be
% used to get the template number in case we want.
basePath='/Users/fred/matlabWork/Yunhui2/VesicleFinderData/120711/';
baseName='004_sq02_1_04';
load([basePath baseName 'rscc.mat']);
% contains:'mxCC','mxVars','mxVesInds','mxTemplInds','m0','m1','partRadius','ds';
return

load([basePath 'Info/' baseName 'mi.mat']);
% mi=rsSortVesicles(mi);  % Sort according to position

% % partRadius=12;
% partRadius=12;
minAmp=0.75;
maxAmp=1.1;
maxVar=40;

n=size(mxCC,1);
nb=2*partRadius+4;  % size of the blanking mask
nb2=nb/2;
blankMask=1-fuzzymask(nb,2,partRadius,1);

results=zeros(1000,5);
ctrs=zeros(1000,2);
vesIndex=single(zeros(1000,1));
k=0;
nFound=0;

disp('Finding particles');
nctr=ceil((n+1/2));
mxCC2=mxCC;  % Search for cc peaks
[amp ix iy]=max2d(mxCC2);
while amp>minAmp
    ccVar=mxVars(ix,iy);
    iv=mxVesInds(ix,iy);
    result=[ iv amp round(ccVar) [ix iy]+round(nctr)];
    k=k+1;
    results(k,:)=result;
    if (amp > minAmp && amp < maxAmp && ccVar < maxVar)
        nFound=nFound+1;
        ctrs(nFound,:)=[ix iy];
        vesIndex(nFound)=iv;
        amps(nFound)=amp;
    end;
    % Blank the vicinity of the found peak
    mxCC2=Mask(mxCC2,[ix iy],blankMask);
    [amp ix iy]=max2d(mxCC2);
    if mod(k,10)==0
        imacs(mxCC2);
        title(amp);
        drawnow;
    end;
end;
results=results(1:nFound,:);
ctrs=ctrs(1:nFound,:);
for i=1:nFound
    mi.particle.x(i)=(ctrs(i,1)-1)*ds;  % should be interpolated!
    mi.particle.y(i)=(ctrs(i,2)-1)*ds;
    mi.particle.vesicle=vesIndex(1:nFound);
end;



%
figure(3);
ShowImageAndBoxes(mxCC,ctrs,20,2,[1 1 0]);

figure(2); clf;
ShowImageAndBoxes(m1,ctrs,20,2,[1 1 0],1e-3);

figure(1);
ShowImageAndBoxes(m0,ctrs,20,2,[1 1 0],1e-3);
