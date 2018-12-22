function [epaVals,ctfVals,ctfImage]=rtReadGctfLogs(mi)
% function [epaVals,ctfVals,ctfImage]=rtReadGctfLogs(mi)
% After Gctf has put various log files into the output directory (mi.imagePath)
% we read the EPA.log, gctf.log and .ctf image and interpret them.
% epaVals is a struct of vectors. ctfVals is a struct of final values that
% Gctf writes in gctf.log. ctfImage is the image in the .ctf file.

% mrcFilename='Aug10_14.58.56_DW.mrc';
mrcFilename=[mi.imagePath mi.imageFilenames{1}];
[pa,nm,ex]=fileparts(mrcFilename);
path=AddSlash(pa);
logNameEPA=[path nm '_EPA.log'];
logNameGctf=[path nm '_gctf.log'];
mrcNameCTF=[path nm '.ctf'];
% logNameGctf='Temp/GctfLog.txt'; %######

ctfImage=zeros(512,512,'single');
if exist(mrcNameCTF,'file')
    [ctfImage,s]=ReadMRC(mrcNameCTF);
end;

fEPA=fopen(logNameEPA);
data=cell(6,1);  % default are null arrays of values.
if fEPA>0 % We got a valid file
    %     nCols=5;  % hard-wired number of columns
    header=textscan(fEPA,'%s',1,'delimiter','\n'); % 1st line of the file
%    disp(header{1});
    data=textscan(fEPA,'%f%f%f%f%f');
    fclose(fEPA);
end;
epaVals=struct;
epaVals.resolution=data{1};
epaVals.ctfSim=data{2};
epaVals.epaRaw=data{3};
epaVals.epaBkgSub=data{4};
epaVals.ccc=data{5};

% by default, return an empty structure for ctfVals.
ctfVals=struct;

fGctf=fopen(logNameGctf);
if fGctf>0 % successful file opening
    gData=textscan(fGctf,'%s','delimiter','\n');
    gLines=gData{1};
    if numel(gLines)>100 % check for approximately full-length file.
        % We grep for the last occurrence of 'LAST CYCLE'
        cPtrs=strfind(gLines,'LAST CYCLE');
        ng=numel(gLines);
        k=ng;
        while k>0 && numel(cPtrs{k})==0
            k=k-1;
        end;
        % Then grep for first occurrence of Defocus_U following LAST CYCLE
        if k>0 % we found the last line containing 'LAST CYCLE'; search subset
            qLines=gLines(k:ng); % cell array of last lines
            qPtrs=strfind(qLines,'Defocus_U');
            nq=numel(qPtrs);
            for j=1:nq
                if numel(qPtrs{j}>0)
                    break;
                end;
            end;
            if j<nq % found the headers and values
                headers=textscan(qLines{j},'%s');
                headers=headers{1};
                vals=textscan(qLines{j+1},'%f');
                vals=vals{1};
                % pick up the Resolution limit
                [headers,vals]=ExtractFinalValue(qLines,'RES_LIMIT',headers,vals);
                [headers,vals]=ExtractFinalValue(qLines,'B_FACTOR',headers,vals);
                ctfVals=struct;
                for i=1:numel(vals)
                    if numel(headers{i})>0 % there is a header name
                        ctfVals.(headers{i})=vals(i);
                    end;
                end;
                
                
            end;            
        end;     
    end;
    fclose(fGctf);
else
    error([logNameGctf ' not found.']);
end;

end


    function [hdrs,dat]=ExtractFinalValue(lines,searchString,hdrs,dat)
        rPtrs=strfind(lines,searchString);
        for i=1:numel(lines)
            r={'';0};
            if numel(rPtrs{i})>0
                r=textscan(lines{i}(rPtrs{i}(1):end),'%s%f');
                break;
            end;
        end;
        hdrs(end+1,1)=r{1};
        dat(end+1,1)=r{2};
        
    end
