function classMembers=mgReadClsFile(filename)
% function classMembers=mgReadClsFile(filename)
% Scan an imagic *.cls file to decode the class members
% per class.  The results are placed into the cell array Cls.
%             classMembers{classnum}= members
% Default filename is classes.cls

classMembers={};

if nargin<1
    filename='classes.cls';
end;

f=fopen(filename,'rt');
if f<=0
    error('File could not be opened');
    return
end;

theLines = textscan(f,'%s','delimiter','\n');
theLines = theLines{1};  % textscan returns a 1x1 cell containing a cell array.
fclose(f);

nl=numel(theLines);
classNum=0;
ind=0;
while ind<nl
    count=0;
    while count==0 && ind<nl
        ind=ind+1;
        line = theLines{ind};
        [m, count]=sscanf(line, '%d');
    end;
    classNum=classNum+1;
    if m(1) ~= classNum
        error(['Should be [classNum numMembers z]: ' line]);
    end;
    nMembers=m(2);
    n=0;
    mem=[];
    while n<nMembers && ind<nl
        ind=ind+1;
        line = theLines{ind};
        [m, count]=sscanf(line, '%d');
        mem=[mem;m];
        n=n+count;
    end;
    classMembers{classNum,1}=mem;
end;
