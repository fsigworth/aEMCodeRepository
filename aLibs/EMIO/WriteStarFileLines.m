function WriteStarFileLines(fields,lines,dataName,fileName)
% function 

% Write a STAR file from the cell arrays fields and lines, which contain
% the fieldname strings of a star file, followed by the data strings.  This
% is the reverse of 
% [fields,lines,dataName]=ReadStarFileLines(name)

nFields=numel(fields);
if nFields<1
    error('Input structure has no fields');
end;
nLines=numel(lines);

%%

fi=fopen(fileName,'w');

fprintf(fi,'\ndata_%s\n',dataName);
fprintf(fi,'\n');

fprintf(fi,'loop_\n');
for i=1:nFields
    fprintf(fi,'_%s\n',fields{i});
end;

for iLine=1:nLines
    fprintf(fi,'%s\n',lines{iLine});
end;
fprintf(fi,'\n');
fclose(fi);

