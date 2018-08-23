function [fields,lines,dataName]=ReadStarFileLines(name)
% Read a Relion .star file and put its contents into two cell arrays. 
% At present we only know how to read a 'loop_' block.
% the fields array contains the strings identifying the data fields for the
% loop.  The lines array contains each line of the loop.
% The dataName is the text following the data_ header line.
%

fields={};
lines={};
dataName='';

fi=fopen(name);
ok=fi>0;

% while ok    
% Skip to the 'data_' line
s=fgetl(fi);
while ~feof(fi) && ~strncmp(s,'data_',5)  % data block
    s=fgetl(fi);
end;

% Better be 'data'
if strncmp(s,'data_',5)  % data block
    dataName=s(6:end);
else
    return
end;

% Check to see if dataName ends with _x where x is numeric,
p=strfind(dataName,'_'); % search for underscores in the remaining name
%     Check that all remaining characters are numeric
index=0;
if numel(p)>0 && p(end)>1
    index=str2num(dataName(p(end)+1:end));
end;
dataName=dataName(1:p(end)-1); % we've stripped the dataName



s=fgetl(fi);
while ~feof(fi) && ~strncmp(s,'loop_',5)
    s=fgetl(fi);
end;

if strncmp(s,'loop_',5)
    s1='_';
    nFields=0;
    fields={};
    ok=1;
    while ok  % pick up fieldnames
        s1=fgetl(fi);
        if s1(1)=='_'
            nFields=nFields+1;
            fields{nFields}=strtok(s1(2:end));
        else
            ok=0;
        end;
    end;
else
    error('No loop statement found');
end;

% at this point s1 contains the first line of data.
lines={};
id=0;
ok=1;
while ok  % pick up data as strings
    if numel(s1)>0 && ~feof(fi)
        id=id+1; % data index (=line number)
        lines{id}=s1;
        s1=fgetl(fi);
    else
        ok=0;
    end;
end;

fclose(fi);
