function [base, final]=ParsePath(pathText)
% function [base, final]=ParsePath(pathText);
% From a path string native to the running machine,
% extract the final directory from it.
% For example, ParsePath('/Users/fred/data') returns
% base='/Users/fred/' and final='data/'.  Both returned strings (if not null
% strings) end with a file separator.

% defaults
pathText=AddSlash(pathText);
    final=pathText;
    base='';
nt=numel(pathText);
p=strfind(pathText,filesep);
np=numel(p);
if p(np)==nt  % ignore trailing slash
    np=np-1;
end;

if np>0
    final=pathText(p(np)+1:nt);
    final=AddSlash(final);
    base=pathText(1:p(np));
end;
