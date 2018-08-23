% ModifyMembraneDensity.m
% Reduce the TM density of the Kv channel to half.
cd('/Users/fred/aEMCodeRepository/AMPAR')
[m,pixA]=ReadEMFile('KvMap.mrc');

m=shiftdim(m,2);

m(:,:,65:100)=m(:,:,65:100)*.5;
% m(:,:,1:45)=0;
ShowSections2(m);

mbnOffsetA=48;
pixA=2;

map=m;

save KvMap.mat map pixA mbnOffsetA
