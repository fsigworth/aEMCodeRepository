function miNames=FindInfoFiles
% For batch processing, this provides the list of file names containing 'mi.'
% found in Info/
% The names (including path) are returned in the cell array miNames.
maxSize=inf; %%%%%% 

d=dir('Info/');
names=cell(numel(d),1);
sizes=zeros(numel(d),1);
for i=1:numel(d)
    names{i}=['Info/' d(i).name];
    sizes(i)=d(i).bytes;
end;
ptr=strfind(names,'mi.');
if isa(ptr,'cell')
    cptr=ptr;
    ptr=zeros(numel(cptr),1);
    for i=1:numel(cptr)
        q=cptr{i};
        if numel(q)>0
            ptr(i,1)=q(end);
        end;
    end;
end;
ok=(ptr>2) & sizes<maxSize;
miNames=names(ok);
