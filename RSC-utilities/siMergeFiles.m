% siMergeFiles.m
% Scan a set of si files and merge them into one huge one.
%

load siNames.mat
nsi=numel(siNames);
allSis=cell(nsi,1);
for i=1:nsi
    disp(siNames{i});
    si1=load(siNames{i});
    allSis{i}=si1.si;
end;

% %% Get all the defocus values and set activeFlags
% 
% doWriteSi=0;
% doResetActiveFlags=1;
% minDefocus=[2 2.5]
% 
% allDefoci=[];
% for i=1:nsi
%     si=allSis{i};
%     if doResetActiveFlags
%         % reset to the original active flags
%         si.activeFlagLog=si.activeFlagLog(1);
%         si.activeFlags=si.activeFlags(1:numel(si.miIndex));
%         si.activeFlags=si.activeFlags(:);
%     end;
%     np=numel(si.miIndex);
%     siDefoci=zeros(np,1);
%     for j=1:np
%         siDefoci(j)=si.mi{si.miIndex(j)}.ctf(1).defocus;
%     end;
%     for k=1:numel(minDefocus)
%         si.activeFlags=[si.activeFlags (siDefoci>=minDefocus(k))];
%         si.activeFlagLog=[si.activeFlagLog; [date '  siScanFiles defocus >= ' num2str(minDefocus(k))]];
%     end;
%     allSis{i}=si;
%     if doWriteSi
%         disp(['Writing ' siNames{i}]);
%         save(siNames{i},'si');
%     end;
%     for k=1:numel(si.activeFlagLog)
%     disp([num2str(sum(si.activeFlags(:,k))) '  ' si.activeFlagLog{k}])
%     end;
%     
%     allDefoci=[allDefoci; siDefoci];
%     disp(' ');
% end;
% numDefoci=numel(allDefoci);
% lims=numDefoci;
% for k=1:numel(minDefocus)
%     lims(k+1,1)=sum(allDefoci>=minDefocus(k));
% end;
% disp('Totals:');
% disp(lims)
% 
% binWidth=0.1;
% figure(1);
% clf;
% hist(allDefoci,binWidth/2:binWidth:max(allDefoci));
% xlabel('Defocus, um');
% ylabel('Frequency');

%% Merge si structures
% --we also add the stackName and stackIndex fields.
%
stackSuffixIn='stack.mrc';
stackSuffixOut='stack.mrcs';
stackUSuffixIn='ustack.mrc';
stackUSuffixOut='ustack.mrcs';

siMerged=allSis{1};

[~,name,~]=fileparts(siNames{1});
siMerged.baseNames={name(1:end-6)}; % delete('tsi.mat')

np0=numel(siMerged.miIndex);
siMerged.mi=siMerged.mi(:);
siMerged.siIndex=ones(np0,1,'int32');
theFields=fieldnames(siMerged);
ignoreFields={'pixA' 'pars' 'mbnOffset' 'weights' 'activeFlagLog' ...
    'ctfs' 'ctf1s' 'pwfs' 'baseNames' 'mergeMode'};
keeps=true(numel(theFields,1));
for i=1:numel(theFields)
    keeps(i) = ~any(strcmp(theFields{i},ignoreFields));
end;
theFields=theFields(keeps);

for i=2:nsi
    si=allSis{i};
    
    si.mi=si.mi(:);
    nmiSum=numel(siMerged.mi);
    si.miIndex=si.miIndex+nmiSum;
    
    np0=numel(si.miIndex);
    si.siIndex=i*ones(np0,1,'int32');
    
    for j=1:numel(theFields)
        siMerged.(theFields{j})=[siMerged.(theFields{j,:}); si.(theFields{j,:})];
    end;
    
    [~,name,~]=fileparts(siNames{i});
    siMerged.baseNames{i,1}=(name(1:end-6)); % delete('tsi.mat')

end;

 