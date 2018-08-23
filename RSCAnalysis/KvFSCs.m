% KvFSCs.m

% cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon128Ld');
cd('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon112ct/mrc')

% [v20t1,s]=ReadMRC('141106v1t1i20.mrc');
% [v20t2,s]=ReadMRC('141106v1t2i20.mrc');
[v20t1,s]=ReadMRC('i30av01.mrc');
[v20t2,s]=ReadMRC('i30bv01.mrc');
[v20tc,s]=ReadMRC('i30cv01.mrc');
% v14t2=ReadMRC('141106v1t2i13.mrc');
n=size(v20t2,1);
ct=n/2+1;

figure(5);
ShowSections2(v20t1);
title('v20t1');
figure(6);
% ShowSections2(v20t2);


% v20t2r=rsRotateImage(v20t2,45*ones(n,1));
v20t2r=reAlignVolumes(v20t1,v20t2);
% v20t2r=v20t2;
ShowSections2(v20t2r);
title('v20t2r');


fsc0=FSCorr(v20t1,v20t2);

%%
bMask=fuzzymask(n,3,18,5,[ct ct 37]);
bMask=fuzzymask(n,3,40,5,[ct ct 60]);
% bMask=1;
figure(7);
ShowSections2(v20t2r.*bMask);
%
fsc1=FSCorr(v20t2r.*bMask,v20t1.*bMask);

f=(1:ct-1)/(n*s.pixA);
figure(1);
plot(f,[fsc0 fsc1 fsc1*0+.143 fsc1*0]);

% WriteMRC(v20t1+v20t2,s.pixA,'141106v1sum.mrc');
