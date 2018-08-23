% rspSortMiList.m
% From stored allMis.mat and allNames.mat, construct the allNamesSorted
% cell array and save it.

infoPath='Info/'; % assumes we're in the experiment directory
% MiLoadAll;  % loads allMis, allNames

% disp('Loading allMis.mat');
% load([infoPath 'allMis.mat']);
% disp('Loading allNames.mat');
% load([infoPath 'allNames.mat']);
% nmi=numel(allMis);
nmi=numel(allNames);

mip=zeros(nmi,1);
defoci=zeros(nmi,1);

for i=1:nmi
    mi=allMis{i};
    if isfield(mi,'ctf') && isfield(mi.ctf,'defocus') && numel(mi.ctf)>1
        defoci(i)=mi.ctf(1).defocus;
    end;
end;
[defVals,ptrs]=sort(defoci);

subplot(211);
plot(defoci);
subplot(212);
plot(defVals);
xlabel('Name index');
ylabel('Defocus, \mum');
grid on;
allNamesSorted=allNames(ptrs);
disp('Writing Info/allNamesSorted.mat.');
save([infoPath 'allNamesSorted.mat'],'allNamesSorted');
