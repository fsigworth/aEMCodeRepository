% SearchMiLogs.m
% Flag mi files according to the log entries.

%txt='rsRefineVesicleFits 2017-11-1';
txts={  'MergeImages'
    'VesicleFinder'
    'rsRefineVesicleFits'
    'meInverseFilterAuto'
    'rsPickingPreprocessor' };

theNames=f2FindInfoFiles;

nmi=numel(theNames);
nt=numel(txts);
outNames={};
outFlags=zeros(nmi,1);
for i=1:nmi
    mi=ReadMiFile(theNames{i});
    for j=1:nt % scan for each text entry
        q=(strfind(mi.log,txts{j}));
        for k=1:numel(q)
            if numel(q{k})>0 && outFlags(i)==j-1 % log{k} has that entry
                outFlags(i)=j;
            end;
        end;
    end;
    %     for j=1:numel(p)
    %         found=found || numel(p{j}>0);
    %     end; if found
    % %        disp(theNames{i});
    %     else
    %         k=k+1; outNames{k}=theNames{i};
    %     end;
    if mod(i,100)==1
        fprintf('.');
    end;
end;
fprintf('\n');
% disp([num2str(k) ' files are missing ''' txt '''.']); disp('The list is
% contained in cell array ''theNames''');
%%
txts{nt+1}='no flag';
disp(' ');
disp('Count  Log entry');
for i=1:nt+1
    j=mod(i,nt+1);
    disp([num2str(sum(outFlags==j)) '  ' txts{i}]);
end;

