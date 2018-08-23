% f2MakeDisplayJpegs.m
    doDeleteFiles=1;
    
    [names, pathName]=uigetfile('*DisDat.mat','Select display data files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(names)
        names={names};
    end;
    cd(pathName);

    %%
    nFiles=numel(names);
    figure(1);
    
    for i=1:nFiles
        name=names{i};
        disp(name);
        q=load(name);
        qNames=fieldnames(q);
        s=q.(qNames{1});  % pick up the first variable in the file
        if ~isa(s,'struct')
            warning([name ' doesn''t contain a struct']);
            break
        end;
        sNames=fieldnames(s);
        nPanels=numel(sNames);
        nr=floor(sqrt(nPanels));
        nc=ceil(nPanels/nr);
        figSize=[250*nc 235*nr];
        disp(figSize);
        pos=get(gcf,'position');
        pos(3:4)=figSize;
        set(gcf,'position',pos);
        
        DrawFigureFromData(s);
        
        p=strfind(name,'DisDat');
        if numel(p)>0
            p=p(1);
            outName=[name(1:p-1) '.jpg'];
            set(gcf,'paperpositionmode','auto');
            print('-djpeg','-r150',outName);  % save the figure.
            if doDeleteFiles
            delete(name);  % delete the original file
            end;
        end;
    end;
