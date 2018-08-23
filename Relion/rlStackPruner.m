% StackPruner
% deselect stack images that are out of bounds

stackSuffix='tstack.mrcs';
activeFlagIndex=inf;

% Get a tsi file to combine
[siName, pa]=uigetfile('*tsi.mat','Select tsi file');
[rootPath, stackPath]=ParsePath(pa);
cd(rootPath);

load([stackPath siName]);
%%
p=strfind(siName,'tsi.mat');
stackName=[stackPath siName(1:p-1) stackSuffix];
disp(['Reading ' stackName]);
stack=ReadMRC(stackName);
disp(' done.');

figure;
clf
hist(stack(:),1000);

%%
[n,n1,nim]=size(stack);
means=mean(reshape(stack,n^2,nim));
sds=std(reshape(stack,n^2,nim));

figure(1);
subplot(211);
plot(sds);
subplot(212);
plot(means);
drawnow;
%
figure(2);
hist(sds,200);
drawnow;
%
sThreshL=median(sds)-0.075;
sThreshU=median(sds)+.075;
mThreshL=-.1;
mThreshU=.1;
bad=sds<sThreshL | sds>sThreshU | means<mThreshL | means>mThreshU;
sum(bad)

pause

%%
afIndex=min(activeFlagIndex,size(si.activeFlags,2)+1)
si.activeFlags(:,afIndex)=~bad;
si.activeFlagLog{afIndex}=[date '  StackPruner'];
sum(si.activeFlags)
si.activeFlagLog

%%
s=input('Write the stack? ','s');
if s=='y'
save([stackPath siName],'si');
disp([stackPath siName ' written.']);
else
    disp('not written.');
end

