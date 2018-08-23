% SetMiParameter
% Set a field in a group of mi files.

field='gainRefName';
value='movies/movie_frames/CountRef_May09_12.59.08.dm4';

% Have the user select some mi files
[fname, pa]=uigetfile('*mi.mat','Select mi files to examine','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;
%%
cd(pa);
nim=numel(fname);
val=zeros(nim,1);
for i=1:nim
    load(fname{i});
    mi.(field)=value;
    save(fname{i},'mi');
    disp(fname{i});
end;
