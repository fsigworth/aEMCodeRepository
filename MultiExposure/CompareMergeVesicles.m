% CompareMergeVesicles

mia=cell(0);
load('/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/Info_100/sq02_10006mi.mat')

mia{1}=mi;
load('/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/Info_010/sq02_10006mi.mat')
mia{2}=mi;
load('/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/Info_110/sq02_10006mi.mat')
mia{3}=mi;
load('/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/Info/sq02_10006mi.mat')
mi.weights=[1 1 1];
mia{4}=mi;
mi.weights=[0 0 1];
mia{5}=mi;

nim=numel(mia);
figure(1);
xs=0:.0001:.01;
for i=1:nim
    subplot(nim,1,i);
    h=hist(mia{i}.vesicle.s/mia{i}.doses(1),xs);
    bar(xs,h);
    ylabel(num2str(mia{i}.weights));
end;
%%
figure(2);
Hs=[];
for i=1:nim
    mi=mia{i};
    mi.ctf(1).B=0;
    mi.ctf(2).B=0;
    mi.ctf(3).B=0;
    mi.doses=[1 1 1];
    H=meGetEffectiveCTF(mi,1024);
    Hs(:,i)=sectr(H);
    labels{i}=num2str(mi.weights);
end;
plot(Hs);
legend(labels);