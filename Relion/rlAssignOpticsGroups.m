% rlAssignOpticsGroups
% Reads a movies.star file (maybe someday a micrographs.star also)
% And assigns optics group numbers based on shot coordinates in the movie
% name, e.g.  .....0053_X+1Y+1-0.tif
% We detect this 6-character string by searching for the last 'X' in the name.
%
inStarName='movies.star';
inStarName='corrected_micrographs.star';
inStarName='micrographs_sub_ctf.star';
[pa,nm,ex]=fileparts(inStarName);
outStarName=[nm '_optics' ex];


[nm,da]=ReadStarFile(inStarName);
d=da{2};

% find the last X in each name
% names=d.rlnMicrographMovieName;
names=d.rlnMicrographName;
startInds=strfind(names,'X');
nl=numel(names);
% pick up the 6-character position strings
posStrings=cell(nl,1);
for i=1:nl
    start=startInds{i}(end);
    % Pick up the 6-character position string
    posStrings{i}=names{i}(start:start+5);
end;
[uStrings,uFirst,uGroups]=unique(posStrings);
ng=numel(uStrings);
% Assign the optics groups
d.rlnOpticsGroup=uGroups;

% Now expand the optics group block, assuming that we have one entry.
op=da{1};
fNames=fieldnames(op);
op1=struct;
% gpNames=cell(ng,1);
% for i=1:ng
%     gpNames{i}=sprintf('OpticsGroup%02d',i);
% end;
gpNames=uStrings;

for i=1:numel(fNames)
    field=fNames{i};
    if strcmp(field,'rlnOpticsGroup')
        op1.(field)=(1:ng)';
    elseif strcmp(field,'rlnOpticsGroupName')
        op1.(field)=gpNames;
    else
        op1.(field)=repmat(op.(field)(1),ng,1);
    end;
end;
da{1}=op1;
da{2}=d;
disp(['Writing ' outStarName]);
WriteStarFile(nm,da,outStarName);

