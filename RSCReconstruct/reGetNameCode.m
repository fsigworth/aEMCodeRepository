function name=reGetNameCode(ri,iter,iTwin,iGroup,version)
%  create the prefix for a shared file, encoded according to group and twin
% like these, for iter=1, iTwin=2 and <bn> is ri.outBasename:
% <bn>i01b  (no group specified)
% <bn>i01b02 (iGroup=2)
% Special case: if iGroup is negative, it's treated as a volume number
% <bn>i01bv02a
if nargin<4
    iGroup=0;
end;
if nargin<5
    version=2;
end;
tChar=char(iTwin+96);
switch version
    case 1
        if iGroup>0
            name=sprintf('%si%02dg%02d%s',ri.outBasename,iter,iGroup,tChar);
        elseif iGroup==0
            name=sprintf('%si%02d%s_',ri.outBasename,iter,tChar);
        else  % negative iGroup means it's actually iVol
            name=sprintf('%si%02dv%02d%s',ri.outBasename,iter,-iGroup,tChar);
        end;
    case 2
        if iGroup>0
            name=sprintf('%si%02d%s%02d',ri.outBasename,iter,tChar,iGroup);
        elseif iGroup==0
            name=sprintf('%si%02d%s_',ri.outBasename,iter,tChar);
        else  % negative iGroup means it's actually iVol
            name=sprintf('%si%02d%sv%02d',ri.outBasename,iter,tChar,-iGroup);
        end;
end;
