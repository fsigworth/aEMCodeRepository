function s=LoadStruct(matFilename)
% function s=LoadStruct(matFilename)
% From a Matlab .mat file, transparently load either a set of variables or
% contents of a struct into fields of the output s.  Thus
% si=LoadStruct('siName.mat') always loads the struct si whether it was
% saved as
% save('siName.mat','si') %  conventional way,
% or:
% save('siName.mat','-struct','si')  % new preferred way to save data.

s=load(matFilename);
fields=fieldnames(s);
if numel(fields)==1 && isstruct(s.(fields{1}))
    s=s.(fields{1});
end;
