% rlMakeKvSolMask.m

cd('/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/Stack2')
[msk,s]=ReadMRC('KvMask192.mrc');

% Double sphere mask for membranes
sph=fuzzymask(192,3,120,6,[96 96 0]) + fuzzymask(192,3,120,6,[96 96 270]);
sphc=1-sph;
SolMsk=min(sphc+msk,1).*fuzzymask(192,3,92,8);
ShowSections(SolMsk,[],45);
WriteMRC(SolMsk,s.pixA,'KvSolMask');



