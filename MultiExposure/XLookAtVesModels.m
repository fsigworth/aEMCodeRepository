% MultiExposure/XLookAtVesModels.m

[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
% pa
if ~iscell(fname)
    fname={fname};
end;
nim=numel(fname);
%%
models=[];

for i=1:nim
    load([pa fname{i}]);
    models(:,i)=mi.vesicleModel;
end;
plot(models)
hold on
plot(mean(models,2),'k-')
hold off;
legend(num2str((1:nim+1)'));