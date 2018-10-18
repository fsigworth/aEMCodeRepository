function [sequence,datenums]=miDecodeLog(mi);
% From the given mi structure, return a vector of 8 indices, whose elements
% correspond to the programs
% 1 Find jump
% 2 Drift tracker
% 3 Merge
% 4 Vesicle [finder]
% 5 RefineVesicleFits
% 6 InverseFilter
% 7 RefineMembraneModel
% 8 PickingPreprocessor
%
% If any of these programs appear in the log, the index will be nonzero.
% Instead it will be the index of the last appearance in the logs.

searchStrings={'k2FindDefocusJump'; 'k2DriftTracker' ; 'MergeImages' ; ...
    'Vesicle' ; 'rsRefineVesicleFits' ; 'meInverseFilter' ; 'RfineMmbraneModelXX' ; ...
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
if sequence(4)<1
% patch for old files where no log entry is made
if isfield(mi,'vesicle') && isfield(mi.vesicle,'x') && numel(mi.vesicle.x)>0
    sequence(4)=sequence(3)+.5;
end;
end;

if nargout>1 % look up timestamps
    nlog=numel(mi.log);
    datenums=zeros(ns,1);
    for i=find(sequence>0)
        strings=split(mi.log{floor(sequence(i))});
        try
            datenums(i)=datenum([strings{end-1} ' ' strings{end}]);
        end;
    end;
end;
