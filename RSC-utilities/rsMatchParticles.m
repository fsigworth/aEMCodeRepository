% rsMatchParticles
% Match found particles with true particle positions.  Make a histogram of
% amplitude values.

outName='amps20m22.mat';
outName='amps30m22.mat';
ampScaler=1;
useTrue=1;

% outName='amps40m1.mat';
% ampScaler=1e3;

[fname, pa]=uigetfile('*mi.mat','Select the mi file to search');
load([pa fname]);
basePath=ParsePath(pa);
%%
distThreshold=15;  % original pixels

% autopicks0=mi0.particle.picks(:,3)==32;  % auto-picked flag
% locs0=mi0.particle.picks(autopicks0,1:2);
%

    autopicks=mi.particle.picks(:,3)==32;
    amps=mi.particle.picks(:,5);

if useTrue && isfield(mi.particle,'true')
    locs0=mi.particle.true(:,1:2);  %%%%
    
    locs=mi.particle.picks(:,1:2);
    nlocs=size(locs,1);
    
    minDist=zeros(nlocs,1);
    for i=1:nlocs
        dists=hypot(locs0(:,1)-locs(i,1),locs0(:,1)-locs(i,1));
        minDist(i)=min(dists);
    end;
    
    ok=minDist<=distThreshold;
    ok=ok & amps>.32;  %%%%
    allowFalse=amps<.34;
    
else
    
    ok=amps>.3;  %%%%%
end;

mi.particle.picks(:,7)=ok;
amps=amps*ampScaler;
xMax=.6;
binWidth=.015;
bins=binWidth/2:binWidth:xMax;
ampsTrue=amps(ok&autopicks);
ampsFalse=amps(~ok & autopicks & allowFalse);

h1=hist(ampsTrue,bins)';

h0=hist(ampsFalse,bins)';
firstBin=find(h0,1);
h0(firstBin)=h0(firstBin+1);
% [binVal, lowBin]=max(h0);
% h0(1:lowBin-1)=h0(lowBin)+h1(lowBin);

bar(bins,[h1 h0],'stacked');
save([basePath outName], 'ampsTrue', 'ampsFalse', 'h0', 'h1', 'bins');
