% rlSelectIOParticles.m
% Given a stack file, set the active flags for inside-out or right-side-out
% particles.
% mi.particle.picks(:,7) is the rso flag.

pickRSOs=false;

% % Put up a file selector for files *si.mat,
%     disp('Getting si file names');
%     [fname, siPath]=uigetfile('*si.mat','Select si file');
%     if isnumeric(siPath)  % user clicked Cancel
%         return
%     end;
%     pa=ParsePath(siPath);  % go to main directory
%     cd(pa);
%     siPath=AddSlash(siPath);
% % Get the si structure
%     load([siPath fname]);
    
np=numel(si.miIndex);
oks=false(np,1);
    for i=1:np
        mi=si.mi{si.miIndex(i)};
        rso=mi.particle.picks(si.miParticle(i),7);
        oks(i)=(rso==pickRSOs);
    end;
    disp([num2str(sum(oks)) ' particles with RSO= ' num2str(pickRSOs) ' out of ' num2str(np)]);
%     naf=size(si.activeFlags,2)+1;
naf=3;
    si.activeFlags(:,naf)=si.activeFlags(:,naf-1) & oks;
    disp(sum(si.activeFlags,1));
    si.activeFlagLog{naf,1}=[date '  rlSelectIOParticles: RSO=' num2str(pickRSOs)];
    disp(si.activeFlagLog);
    
    
%%    
%     Trim mi particles, deleting all the background flags

nmi=numel(si.mi);
maxPI=zeros(nmi,1);
numPI=zeros(nmi,1);
for i=1:np
        imi=si.miIndex(i);
        maxPI(imi)=max(maxPI(imi),si.miParticle(i));
end;
figure(3);
plot(maxPI);

for j=1:nmi
    numPI(j)=size(si.mi{j}.particle.picks,1);
    si.mi{j}.particle.picks(maxPI(j)+1:end,:)=[];
end;
plot([maxPI numPI]);