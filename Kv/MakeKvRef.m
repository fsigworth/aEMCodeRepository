% MakeKvRef.m
% Make a scaled reference for Kv reference or mask

refMode=0;

if refMode
    nm='KvRef';
else
    nm='KvMask';
end;

workingDir='/gpfs/ysm/home/fjs2/data/';
outDir=AddSlash(pwd); % current directory
[m0,s]=ReadMRC([workingDir nm '1.05p256.mrc']);

pixA=MyInput('pixel size, A',s.pixA);
n=MyInput('image size',size(m0,1));
m1=DownsampleGeneral(m0,n,s.pixA/pixA);
if ~refMode
    m1=max(0,min(1,m1));
end;
outName=[outDir nm num2str(pixA,4) 'p' num2str(n) '.mrc'];
WriteMRC(m1,pixA,outName);
disp(['wrote ' outName '.' ]);

