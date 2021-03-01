% rlMakeKvSolMask.m

% cd('/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/Stack2')
[msk,s]=ReadMRC('KvMask192.mrc');

% Double sphere mask for membranes
sph=fuzzymask(192,3,120,6,[96 96 0]) + fuzzymask(192,3,120,6,[96 96 270]);
sphc=1-sph;
solMsk=min(sphc+msk,1).*fuzzymask(192,3,92,8);
ShowSections(solMsk,[],45);
n=size(solMsk);
solMsk1=reshape(max(0,min(solMsk(:),1)),n);
% WriteMRC(solMsk1,s.pixA,'KvSolMask');
solMsk21=DownsampleGeneral(solMsk,128,s.pixA/2.1);
n=size(solMsk21);
solMsk2=reshape(max(0,min(solMsk21(:),1)),n);
% WriteMRC(solMsk2,2.1,'KvSolMask2.1A128.mrc');


