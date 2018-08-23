function [names basename defoci ind]=FindHidekiFileSetSeq(d,ind)
% Find defocus pairs just as successive dm3 files.  Odd are close to focus,
% even are far from focus.
ind=ind+1;
defoci=[1 5];
ext='.dm3';
basename='';
names={};
[ind name1]=GetNextEMFile(d,ind,ext);
if ind>0
    [ind name2]=GetNextEMFile(d,ind+1,ext);
end;
if ind>0
    names{1}=name1;
 names{2}=name2;
end;
