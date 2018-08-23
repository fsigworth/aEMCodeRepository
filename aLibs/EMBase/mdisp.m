function mdisp(handles,name,datum)
% function mdisp(handles,datum)
% function mdisp(handles,datum)
% function mdisp(datum)
% Replacement for disp() that prints to one or more files using the mprintf
% function.  handles is an array of file handles.  It calls the function
% mprintf.  A handle of 1 or missing handle argument means print to the
% command window. If name is given, prints that string and then the datum
% value.
% Currently works with strings, numeric up to dim=2,
% structs and struct arrays
if isa(handles,'struct')
    idString=handles.idString;
    handles=handles.handles;
else
    idString='';
end;

nid=numel(idString);
if nid>0
    mprintf(handles,'%s ',idString);
    xspaces=nid+1;
else
    xspaces=0;
end;

switch nargin
    case 1
        datum=handles;
        handles=1;
    case 2
        datum=name;
    case 3
        mprintf(handles,'%s =',name)
        xspaces=numel(name)+2;
end;

sz=size(datum);

if isa(datum,'char')
    mprintf(handles,' %s\n',datum);
elseif isa(datum,'numeric') && sz(1)>0
     mprintf(handles,'   %s\n',num2str(datum(1,:)));
    for irow=2:sz(1)
        mprintf(handles,'%s   %s\n',blanks(xspaces),num2str(datum(irow,:)));
    end;
elseif isa(datum,'struct') && prod(sz)>0
    if numel(datum)>1  % struct array, give values for first element
        sz=size(datum);
        mprintf(handles,'[ %s ] struct array.  First element:\n',num2str(sz));
        datum=datum(1);
        xspaces=0;
    end;
    fields=fieldnames(datum);
    len=0;
    for i=1:numel(fields)
        len=max(len,numel(fields{i}));
    end;
    for i=1:numel(fields)
        nsp=(i>1)*xspaces;
    format=['    %' num2str(len+nsp) 's: %s\n'];
        mprintf(handles,format,fields{i},num2str(datum.(fields{i})));
    end;
else        % all other classes: just give the class type.
    cls=class(datum);
    sz=size(datum);
    mprintf(handles,'[ %s ] %s array\n',num2str(sz),cls);
end;
