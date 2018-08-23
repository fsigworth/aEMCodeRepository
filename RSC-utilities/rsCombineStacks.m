% rsCombineStacks
% Combine stacks from 151117 and 140625

cd('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/StackMerge')

s1Name='sq10_350fp128t';  % from 140625
af1=4;  % active flags
s2Name='sq03_1fp128t';
af2=[5 6];
si1=load([s1Name 'si.mat']);
si1=si1.si;

m1=ReadMRC([s1Name 'stack.mrc']);
m1u=ReadMRC([s1Name 'ustack.mrc']);
n=size(m1,1);
n1=size(m1,3);
%%  % scale this old image set by dose
dose1=zeros(n1,1);
for i=1:n1
    dose1(i)=si1.mi{si1.miIndex(i)}.doses(1);
end;
active1=si1.activeFlags(:,4);  % classes 3 and 4
m1s=m1.*shiftdim(repmat(dose1,1,n,n),1);
m1us=m1u.*shiftdim(repmat(dose1,1,n,n),1);
%% Check that our scaling is good.
disp('Computing spectrum 1');
% sp1=RimSpectrum(m1,msk);
sp1=RadialPowerSpectrum(m1s,1,2);

disp('done.');

%%
fMin=round(n*.06);
fMax=round(n*.2);
var1=mean(sp1(fMin:fMax,:),1)';

plot(var1);

v1=mean(var1(active1))

%%
si2=load([s2Name 'si.mat']);
si2=si2.si;

m2=ReadMRC([s2Name 'stack.mrc']);
m2u=ReadMRC([s2Name 'ustack.mrc']);
n2=size(m2,3);
%%  % scale this old image set by dose.  Due to a bug, we have to scale it
% by dose^1.5
dose2=zeros(n2,1);
for i=1:n2
    dose2(i)=si2.mi{si2.miIndex(i)}.doses(1);
end;
active2=any(si2.activeFlags(:,af2),2);  % classes 4 and 5
m2s=m2.*shiftdim(repmat(dose2.^1.5,1,n,n),1);
m2us=m2u.*shiftdim(repmat(dose2.^1.5,1,n,n),1);
%% Check that our scaling is good.
disp('Computing spectrum 2');
% sp1=RimSpectrum(m1,msk);
sp2=RadialPowerSpectrum(m2s,1,2);
disp('done.');
%

var2=mean(sp2(fMin:fMax,:),1)';

%%
plot(var2);

v2=mean(var2(active2))

%% rescale the second stack to match the first.
scl=sqrt(v1/v2);
m2sc=m2s.*scl;
m2usc=m2us.*scl;


%% Extract the good images

[si1g,m1sg,m1usg]=rsStackSplit(active1,si1,m1s,m1us);
[si2g,m2scg,m2uscg]=rsStackSplit(active2,si2,m2sc,m2usc);

[si3,m3]=rsStackConcatenate(si2g,m2scg,si1g,m1sg);
[si3,m3u]=rsStackConcatenate(si2g,m2uscg,si1g,m1usg);
%%
figure(1);
ImagicDisplay3(BinImage(Crop(m3,96,1),4))

%%
outName='Merge1415t';

si=si3;
disp(['Writing ' outName 'si.mat and stack files...']);
save([outName 'si.mat'],'si');
WriteMRC(m3,si.pixA,[outName 'stack.mrc']);
WriteMRC(m3u,si.pixA,[outName 'ustack.mrc']);
disp('done');



na1=sum(active1);
xStack=m1s(:,:,active1);
na2=sum(active2);
xStack(:,:,na1+1:na1+na2)=m2sc(:,:,active2);



