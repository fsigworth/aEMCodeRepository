% function si=rlSetActiveFlags(si,starName)
% 
nNew=size(si.activeFlags,2)+1;
[nm,pa]=uigetfile('*.star');
starName=[AddSlash(pa) nm];
disp(['Reading ' starName]);

[~,dat,ok]=ReadStarFile(starName);
if ~ok
    error(['Couldn''t read ' starName]);
end;
d=dat{1};
si.activeFlags(:,nNew)=false;

ni=numel(d.rlnImageName);
for i=1:ni
    ind=rlDecodeImageName(d.rlnImageName{i});
    si.activeFlags(ind,nNew)=true;
end;
si.activeFlagLog{nNew,1}=['rlSetActiveFlags ' starName];
% Show the results
si.activeFlagLog
sum(si.activeFlags)
