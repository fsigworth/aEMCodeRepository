function [val, ok]=rspGetGValue(defaultVal)
% Get a numeric value from keystrokes using ginput().  This does not echo
% to the command line.
ok=0;
str=char;  % a zero-element string
[x y b]=Myginput(1,'circle');
while numel(b)>0 && b>3  % Return or click exits.
%     Test if the character is in the set {+ - . 0..9 e E}
    if any(b==[43 45 46 48:57 69 101])
        str=[str char(b)];
    end;
    [x y b]=Myginput(1,'circle');
end;
newVal=str2double(str);
if numel(newVal)>0
    val=newVal;
    ok=1;
else
    val=defaultVal;
end;