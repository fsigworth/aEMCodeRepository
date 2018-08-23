function [s,ok]=reCheckAndLoadMat(path,fileName,timeout,cycleTime,numFields,plainLoad)
% For reReconstruct, load the mat file and assign fields to the struct s.
if nargin<3
    timeout=0;
end;
if nargin<4
    cycleTime=1;
end;
if nargin<5
    numFields=1;  % minimum number of fields that should be present
end;
if nargin<6
    plainLoad=0;
end;

path=AddSlash(path);

ok=false;
s=[];
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
        try
            nTrials=10;
            complete=false;
            while ~complete && nTrials>0
                if plainLoad
                    s=load([path fileName]);
                else
                    s=LoadStruct([path fileName]);
                end;
                complete=numel(fieldnames(s))>=numFields;
                if ~complete
                    disp([' retry ' num2str(nTrials)]);
                    pause(5);
                end;
                nTrials=nTrials-1;
            end;  % while
            return
        catch
            ok=false;
        end;
    end;
    if ok
        return
    end;
    pause(cycleTime)
    t=toc;
end;
