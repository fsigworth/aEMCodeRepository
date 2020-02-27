function [epaVals,ok]=rtReadEPALog(logName)
epaVals=struct;
fEPA=fopen(logName);
data=cell(6,1);  % default are null arrays of values.
if fEPA>0 % We got a valid file
    %     nCols=5;  % hard-wired number of columns
    header=textscan(fEPA,'%s',1,'delimiter','\n', ...
        'TreatAsEmpty',{'-nan','nan'}); % 1st line of the file
%    disp(header{1});
%     data=textscan(fEPA,'%f%f%f%f%f');
    data=textscan(fEPA,'%f%f%f%f%f%f');
    fclose(fEPA);
    ok=1;
else
    ok=0;
    return
end;
epaVals.resolution=data{1};
epaVals.ctfSim=data{2};
epaVals.epaRaw=data{3};
epaVals.epaBkgSub=data{4};
% epaVals.epaRawMinusBg=data{5};
epaVals.ccc=data{5};
