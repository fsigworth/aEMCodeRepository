% rsDeleteParticles
% Delete all the particle picks from mi files.
doDelete=1;

    disp('Selecting mi files');
    [fname, pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
    cd(rootPath);
    if doDelete
        str='deleted.'
    else
        str='';
    end;
    %%
for i=1:numel(fname)
    name=[infoPath fname{i}];
    mi=ReadMiFile(name);
    nPicks=size(mi.particle.picks,1);
    disp([name ' ' num2str(nPicks) ' ' str]);    
    if doDelete
        mi.particle=struct;
mi.particle.picks=[];
mi.particle.autopickPars=[];
WriteMiFile(mi,name);
    end;
end;

