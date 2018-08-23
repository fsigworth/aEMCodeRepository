% rsPruneHighBobs.m
% Read an si file, check for particles that are more than maxShift beyond
% the expected outside-vesicle position, and create a new set of
% si.activeFlags.

maxShift=6;

[siName,pa]=uigetfile('*si.mat','Select si file to read');
load([AddSlash(pa) siName]);
mis=si.mi;
nFlagSets=size(si.activeFlags,2);
nLogs=size(si.activeFlagLog,2);
for i=1:nFlagSets
    disp([num2str(sum(si.activeFlags(:,i))) '  ' si.activeFlagLog{i}]);
end;
iRefFlags=MyInput('Reference set',nFlagSets);

% Look at the geometry of YClick vs rVesicle+-mbnOffset
si1=si;
np=numel(si1.miIndex);
rsos=false(np,1);
for i=1:np
    rsos(i)=si1.mi{si1.miIndex(i)}.particle.picks(si1.miParticle(i),7);
end;


%%
isoOffset=0;
mbnOffset=si1.mbnOffset-0;  % about 20
sel=rsos;
figure(1);
mo=-mbnOffset;
labels={'RSO particles' 'ISO particles'};
for ipol=1:2
    subplot(2,2,ipol);
    rv=si1.rVesicle(sel);
    yc=si1.yClick(sel);
    plot(rv,yc,'b.',rv,rv+mo,'k-',rv,rv+abs(mo)+maxShift,'r-',rv,rv-mo,'g-');
    title(labels{ipol});
    xlabel('Vesicle radius, pixels');
    ylabel('YClick');
    
    bob=yc-mo-rv;
    sinb=(yc-mo)./rv;  % sin beta
    betas=asin(sinb)*180/pi;
    subplot(2,2,ipol+2)
%     hist(real(betas),100);
%     hist(sinb,100);
    hist(bob,100);
    xlabel('Bob, pixels');
    ylabel('Frequency');
    sel=~sel;
    mo=mbnOffset-isoOffset;
end;

thresh=abs(mbnOffset)+maxShift;
inliers=(si1.yClick-si1.rVesicle<thresh);
numOutOfBounds=[sum(~inliers) sum(si.activeFlags(:,iRefFlags)&~inliers)]

b=input('Add active flags? ', 's');
if b=='y'
    si.activeFlags(:,nFlagSets+1)=inliers & si.activeFlags(:,iRefFlags);
    si.activeFlagLog{nFlagSets+1}=[date '  rsPruneHighBobs'];
    disp(si.activeFlagLog(:));
    disp(['Saving ' siName]);
    save([AddSlash(pa) siName],'si');
end;


