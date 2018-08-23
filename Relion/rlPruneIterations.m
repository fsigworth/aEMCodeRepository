% rlPruneIterations
% % From relion Class2D, Class3D, Refine3D etc. directories, delete the files of all but
% % the last or last two iterations.

doDelete=1;
extraIter=1;  % leave the penultimate iteration

disp('Selecting a Class or Refine directory.');
dirName=uigetdir;
d=dir(dirName);
fNames=cell(0,1);
fIters=zeros(0,1);
fRuns=cell(0,1);

% Extract the run string and the iteration number from each
for i=1:numel(d)
    fNames{i}=d(i).name;
    fRuns{i}='';
    fIters(i)=-1;
    usPtrs=strfind(d(i).name,'_');
    itPtr=strfind(d(i).name,'_it');
    if numel(usPtrs)>0 && usPtrs(1)>1 && numel(itPtr)>0
        fRuns{i}=d(i).name(1:usPtrs(1)-1);
        fIters(i)=sscanf(d(i).name(itPtr(1)+3:end),'%d');
    end;
end;

%%

runString='';
j=1;
k=0;
rPtrs=[];
while j<numel(d)
    while j<numel(d) && strcmp(fRuns{j},runString)
        j=j+1;
    end;
    runString=fRuns{j};
    k=k+1;
    rPtrs(k)=j;
    if j>1 && numel(fRuns{j-1})>0 && k>1 % end of a nontrivial run
        pStart=rPtrs(k-1);
        maxIter=fIters(j-1);
        for ip=pStart:rPtrs(k)-1
            if fIters(ip)<maxIter-extraIter
                delName=fNames{ip};
                str=['rm ' AddSlash(dirName) delName];
                disp(str);
                if doDelete
                    system(str);
                end;
                
            end;
        end;
    end;
        
end;

