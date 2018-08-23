% function rsGctfProcessor(miNames)
load allNames
miNames=allNames;
nmi=numel(miNames);
disp('rsGctfProcessor');
% nmi
% nmi=min(nmi,420) %%%%%
for i=1:nmi
    miName=miNames{i};
    disp(miName);
    mi=ReadMiFile(miName);
    mi=rtGctfRunner(mi);
    WriteMiFile(mi,miName);
end;

    