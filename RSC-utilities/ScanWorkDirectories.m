% ScanWorkDirectories

% load workDirs
nw=numel(workDirs);
checkDirs={'Info/';'Micrograph/';'Merged/';'Jpeg'};
multiplicity=[2; 2; 5; 4];
for i=1:nw
    disp(workDirs{i});
    cd(workDirs{i});
    for j=1:numel(checkDirs);
        d=dir(checkDirs{j});
        nf=0;
        for k=1:numel(d);
           if ~d(k).isdir
               nf=nf+1;
           end;
        end;
        disp([checkDirs{j} ' ' num2str([nf nf/multiplicity(j)])]);
    end;
    disp(' ');
end;
