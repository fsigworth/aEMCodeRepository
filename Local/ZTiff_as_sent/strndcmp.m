function ok=strndcmp(str1,pattern,n)
% function ok=strndcmp(str1,pattern,n)
% Like strncmp, but compares the *last* n characters of the two strings.\
% The argument n is optional; by default it is the length of the pattern.
% returns true if n=0; false if str1 or pattern are shorter than n.

if nargin<3
    n=numel(pattern);
end;
if numel(str1)<n || numel(pattern)<n
    ok=false;
else
    ok=strcmp(str1(end-n+1:end),pattern(end-n+1:end));
end;
