% SearchMiLogs.m
% Find mi files lacking a entry in the log field.

%txt='rsRefineVesicleFits 2017-11-1';
txt='VesicleFinder';

theNames=f2FindInfoFiles;

nmi=numel(theNames);
outNames={};
k=0;
for i=1:nmi
    mi=ReadMiFile(theNames{i});
    p=strfind(mi.log,txt);
    found=0;
    for j=1:numel(p)
        found=found || numel(p{j}>0);
    end;
    if found
%        disp(theNames{i});
    else
        k=k+1;
        outNames{k}=theNames{i};
    end;
    if mod(i,100)==1
        fprintf('.');
    end;
end;
fprintf('\n');
disp([num2str(k) ' files are missing ''' txt '''.']);
disp('The list is contained in theNames');

