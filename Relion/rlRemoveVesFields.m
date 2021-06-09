% rlRemoveVesFields.m
% remove vesicleR, vesiclePsi etc. from any star files in the current
% directory, and write out the modified files as noVes_*

d=dir;
for i=3:numel(d)
    [pa,nm,ex]=fileparts(d(i).name);
    if ~d(i).isdir && strcmp(ex,'.star')
        name=d(i).name;
        disp(['Reading ' name]);
        [nms,das]=ReadStarFile(name);
%%
        ds=das{2};
        fnames=fieldnames(ds);
        inds=strncmp('ves',fnames,3);
        if any(inds)
            dat=das;
            disp('Fields removed:');
            for i=find(inds)'
                disp(fnames{i});
                ds=rmfield(ds,fnames{i});
            end;
            dat{2}=ds;
            outName=['noVes_' name];
            disp(['Writing ' outName]);
            WriteStarFile(nms,dat,outName);
        else
           disp(' no ves fields found.');
        end;
    end;
end;    
    
            