% miCountParticles

% MiLoadAll

nmi=numel(allMis);
isoParticles=0;
rsoParticles=0;

for i=1:nmi
    q=allMis{i}.particle;
    if isfield(q,'picks') && size(q.picks,1)>0
        rsoParticles=rsoParticles+sum(q.picks(:,7)==1);
        isoParticles=isoParticles+sum(q.picks(:,7)==0);
    end;
end;

isoParticles
rsoParticles
