function sequence=miDecodeLog(mi);
% Read an mi file and return a indices, whose elements have the meaning
% 1 Find jump
% 2 Drift tracker
% 3 Merge
% 4 FindVesicles
% 5 RefineVesicleFits
% 6 InverseFilter
% 7 RefineMembraneModel
% 8 PickingPreprocessor

% Then the log index containing last appearance of each type is put in the sequence array.
searchStrings={'k2FindDefocusJump'; 'k2DriftTracker' ; 'MergeImages' ; ...
    '**FindVesicles**' ; 'rsRefineVesicleFits' ; 'meInverseFilter' ; 'RfineMmbraneModelXX' ; ...
    'rsPickingPreprocessor'};
ns=numel(searchStrings);
mask=false(1,ns);
sequence=zeros(1,ns);
for i=1:ns
    str=searchStrings{i};
    nl=numel(str);
    bools=strncmp(mi.log,str,nl);
    if any(bools)
        sequence(i)=find(bools,1,'last');
        mask(i)=true;
    end;
end;
% special case for FindVesicles
if sequence(4)>0
    return
end;
if isfield(mi,'vesicle') && isfield(mi.vesicle,'x') && numel(mi.vesicle.x)>0
    sequence(4)=sequence(3)+.5;
end;

