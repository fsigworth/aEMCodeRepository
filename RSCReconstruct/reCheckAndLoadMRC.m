function [m,s,ok]=reCheckAndLoadMRC(path,fileName,timeout,cycleTime)
if nargin<2
    timeout=0;
end;
if nargin<3
    cycleTime=1;
end;

path=AddSlash(path);
ok=false;
s=struct;
s.err=1;
m=[];
t=-.001;  % go through at least once
tic;
while t<timeout && ~ok
    d=dir(path);
    i=1;
    while i<=numel(d) && ~ok
        ok=strcmp(fileName,d(i).name);
        i=i+1;
    end;
    if ok
%         disp(['Reading ' path fileName]);
        [m,s]=ReadMRC([path fileName]);
        if ~s.err
            ok=true;
            return
        end;
    end;
    pause(cycleTime)
    t=toc;
end;
