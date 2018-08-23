function [iout name]=GetNextEMFile(d,i,ext)
% function [iout name]=GetNextEMFile(d,i,ext)
% find a DM3 file.  Starting with directory d entry i, search for the next 
% file with the given extension e.g. '.dm3' (case insensitive, but includes
% the dot). The returned value iout is the index of the
% found directory entry (or -1 if nothing found).
% The returned value name is the same as d(iout).name.
% To search for consecutive files, be sure to increment the
% index variable between calls to this function.
% Here is an example to get two DM3 files
%     d=dir;
%     [i1 name1]=GetNextEMFile(d,1,'.dm3');
%     [i2 name2]=GetNextEMFile(d,i1+1,'.dm3');
%     if i1>0 && i2>0
%         ... use the files ...
%     end;

iout=i;
name='';
ok=0;
if i<=0
    iout=-1;
    return
end;
for i=i:numel(d)
    name=d(i).name;
    [pa nm ex]=fileparts(name);
    if strcmpi(ex,ext)
        ok=1;
        break;
    end;
end;
if ok
        iout=i;
    else
        iout=-1;
    end;
end
