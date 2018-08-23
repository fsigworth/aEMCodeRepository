function [names basename defoci ind]=FindHidekiFileSet(d,ind, seqMode,numExposures)
% function [names basename defoci ind]=FindHidekiFileSet(d,ind, seqMode,numExposures)

if seqMode==1
    % Here ind is the numeric value of the filename, and starts with 1.
    defoci=[];
    ext='.dm3';
    basename='';
    names={};
%     Find the file whose name starts with ind
    for i=1:numel(d)  % search for the starting file
        name=d(i).name;
        [pa nm ex]=fileparts(name);
        if strcmp(ex, ext)  % Has the right extension
        val=sscanf(nm,'%d');
        if numel(val)>0 && val==ind % Found a numeric name
            names{1}=name;
            basename=name;
            break
        end;
    end;
    done=0;
    for k=1:numExposures
       while ~done && i<numel(d)
           i=i+1;
           name=d(i).name;
           [pa nm ex]=fileparts(name);
        val=sscanf(nm,'%d');
           if strcmp(ex, ext) && numel(val)>0 && val==ind+k
               names{k}=name;
               done=1;
           end;
       end;
    end;

else  % names of the form 001u1000.dm3 and 001u4000.dm3 
    formatString='%3d%1s%d';
    ext='.dm3';
    maxNFiles=3;
    nFound=0;
    names={};
    defoci=[];
    [ind name]=GetNextEMFile(d,ind,ext);
    letters=isletter(name);
    p=find(letters,1,'first');
    basename=name(1:p-1);
    if ind>0
        names{1}=name;
        a1=sscanf(name,'%d%1s%d');  % interpret the leading decimal number
        if numel(a1)>=3
            str=char(a1(2));
            defoci=a1(3);   % defocus in nm
            nFound=1;
            for nFound=1:maxNFiles
                [i2 name]=GetNextEMFile(d,ind+1,ext);
                a2=sscanf(name,'%d%1s%d');
                %         name2=name
                if i2<=0 || numel(a2)<3 || a2(1)~=a1(1)
                    break
                end;
                names{nFound+1}=name;
                defoci(nFound+1)=a2(3);
                ind=i2;  % search for next one
            end;
        end;
    end;
    defoci=defoci/1000;  % convert to um
    % names
    % defoci
end;