function WriteStructText(mi,filename,exceptions)
% Write the structure mi to a text file.  This is used by WriteMiText for
% storing the mi structure for example.

if nargin<2 || numel(filename)<1
    filename='';
    f=1; % write to command window
else
    f=fopen(filename,'w');
end;
if nargin<3
    exceptions.fieldNames={};
    exceptions.rowSize=[];
    exceptions.writeClass={};
end;
% exceptions: exceptions for certain field names.
% rowSize: max number of values per row.  =0 for 'transpose' display
% writeClass: special for double and hex.
% For example:
%     ex.fieldNames={'identifier' 'mergeMatrix' 'picks' 'frameSets' 'data'};
%     ex.rowSize=[       1             9            0         0        60 ]
%     ex.writeClass={'double'          ''           ''        ''     'hex'}

defaultRowSize=8;
maxNumbersInLine=10;
neForRowIndex=100;  % number of elements in an array to turn on row indices.
iquant=4;  % quantum of indentation

% Write a comment line containing the file name and date
[pa,nm]=fileparts(filename);
mprintf(f,'%s\n',['% ' nm ' ' TimeStamp]);  % first line is a comment

% Write each field of mi
fields=fieldnames(mi);
for fIndex=1:numel(fields)      % wrote each field of the struct
    WriteField(mi.(fields{fIndex}),fields{fIndex},0);
end;

if f>1 % we opened a file
    fclose(f);
end;

    function WriteField(val,fieldName,indent)
        c=class(val);
        iEx=find(strcmp(fieldName,exceptions.fieldNames),1,'first');
        ex=struct;
        if numel(iEx)>0
            ex.rowSize=exceptions.rowSize(iEx);
            ex.writeClass=exceptions.writeClass{iEx};
        else
            ex.rowSize=defaultRowSize;
            ex.writeClass='';
        end;
        switch c
            case 'struct'
                WriteStruct(val,fieldName,0,indent);
            case {'single' 'double' 'int16' 'int32' 'uint8' 'logical'}
                ne=numel(val);
                sz=size(val);
                isComplex=any(imag(val(:)));
                if strcmp(ex.writeClass,'double')
                    precision=16;
                    suffix='d';
                else
                    precision=7;  % single precision digits
                    suffix='';
                end;
                if ne<maxNumbersInLine && sz(1)<2 % a single row
                    %                     Turn off the silly warning in v2014b about precision
                    %                     being ignored if val is integer.
                    warning('off','MATLAB:num2str:integerSecondArgument');
                    if ne<1
                        str=[blanks(indent) fieldName ' = []\n'];
                    elseif ne==1
                        str=[blanks(indent) fieldName ' = ' num2str(val(:)',precision) suffix '\n'];
                    else
                        str=[blanks(indent) fieldName ' = [ ' num2str(val(:)',precision) ' ]\n'];
                    end;
                    warning('on','MATLAB:num2str:integerSecondArgument');
                    
                    mprintf(f,str);
                    return
                else
                    if strcmp(ex.writeClass,'hex')
                        str=[blanks(indent) fieldName ' ' SizeString(sz) ' hex:\n'];
                        mprintf(f,str);
                        WriteHex(val,ex.rowSize);  % don't indent
                    else
                        if ex.rowSize>0  % usual case
                            rowSize=ex.rowSize;
                            transChar='';
                        else
                            rowSize=sz(2);
                            val=val';
                            transChar='''';  % single tic
                        end;
                        str=[blanks(indent) fieldName ' ' SizeString(sz) transChar ':\n'];
                        mprintf(f,str);
                        
                        showInd=isComplex || (ne>=neForRowIndex);
                        istr='';
                        ind=0;
                        ptr=1;
                        while ptr<=ne
                            ptr1=min(ne,ptr+rowSize-1);
                               str=sprintf('%12.6g ',real(val(ptr:ptr1)));
                            if showInd
                                ind=ind+1;
                                istr=sprintf('%4d: ',ind);
                            end;
                            mprintf(f,'%s%s%s\n',blanks(indent+iquant),istr,str);
                            if isComplex % print a second line below it.
                               str=sprintf('%12.6g ',imag(val(ptr:ptr1)));
                                istr='   +i ';                                
                                mprintf(f,'%s%s%s\n',blanks(indent+iquant),istr,str);
                            end;
                            ptr=ptr1+1;
                        end;
                    end;
                end;
            case 'char'
                mprintf(f,'%s\n',[blanks(indent) fieldName ' = ' val]);
            case 'cell'
                ne=numel(val);
                sz=size(val);
                if ne>0
                    cls=[' ' class(val{1})];
                else
                    cls='';
                end;
                %                 if ne<maxStringsInLine
                %                     mprintf(f,[blanks(indent) fieldName ' = ' cellstr(val) '\n']);
                %                 else
                mprintf(f,[blanks(indent) fieldName ' ' SizeString(sz) ' cell' cls ':\n']);
                for irow=1:ne
                    if isnumeric(val{irow})
                        nameStr=[fieldName '{' num2str(irow) '}'];
                        WriteField(val{irow},nameStr,indent);
                    else  % a string
                        str=[blanks(indent+iquant) char(val{irow})];
                        mprintf(f,'%s\n',str);
                    end;
                end;
                %                 end;
        end
    end

    function WriteStruct(val,name,index, indent)
        %         index=0 if this is an ordinary call to WriteStruct;
        %         index takes values 1,2,... for subsequent elements of struct
        %         array
        sz=size(val);
        ne=numel(val);
        fnames=fieldnames(val);
        nf2=numel(fnames);
        if index==0
            str=[blanks(indent) name ' ' SizeString(sz,nf2) ' struct:'];
        else
            str=[blanks(indent) name '(' num2str(index) ')'];
        end;
        mprintf(f,'%s\n',str);
        if ne>1  % array of stucts: call recursively
            for i2=1:ne
                WriteStruct(val(i2),name,i2,indent);
            end;
        else  % ordinary struct
            if nf2<1
                return
            end;
            %             See if we have a lot of rows in 2D arrays, with
            %                           multiple arrays having the same numbr of rows.
            numRows=zeros(nf2,1);
            numDims=zeros(nf2,1);
            for i2=1:nf2
                numRows(i2)=size(val.(fnames{i2}),1);
                numDims(i2)=ndims(val.(fnames{i2}));
            end;
            mxRows=max(numRows);
            longRows=(numRows==mxRows & numDims<3);
            if mxRows>2 && sum(longRows)>2  % worth doing columns
                WriteFieldColumns(val,fnames(longRows),indent+iquant);
            else
                longRows=false(1,nf2);
            end;
            fnames=fnames(~longRows); % write the remaining rows.
            nf2=numel(fnames);
            for i2=1:nf2
                WriteField(val.(fnames{i2}),fnames{i2},indent+iquant);
            end;
        end;
    end

    function WriteFieldColumns(val,names,indent)
        %         Write the given fields of the struct val as columns
        indWidth=5;  % size of index column
        indFmt='%4d: ';
        nf2=numel(names);
        nRows=size(val.(names{1}),1);
        nCols=zeros(nf2,1);
        widths=zeros(nf2,1);
        fmt=cell(nf2,1);
        isComplex=false(nf2,1);
        hdrString=[blanks(indent+iquant+indWidth) '>>  ' blanks(indWidth)];
        for ic=1:nf2
            field=val.(names{ic});
            isComplex(ic)=any(imag(field(:))~=0);
            nCols(ic)=size(field,2);
            switch class(field)
                case {'single' 'double' }
                    widths(ic)=13;
                    fmt{ic}='%12.6g ';
                case {'logical'}
                    widths(ic)=2;
                    fmt{ic}='%2d ';
                case ('int16')
                    widths(ic)=4;
                    fmt(ic)='%4d ';
            end;
            str=sprintf(['%-' num2str(widths(ic)*nCols(ic)) 's '],...
                [names{ic} ' ' SizeString(size(field))]);
            hdrString=[hdrString str];
        end;
        mprintf(f,'%s\n',hdrString);
        
        %         Write the rows
        for k=1:nRows
            rowStr=[blanks(indent+iquant+1) sprintf(indFmt,k)];
            for ic=1:nf2
                str= sprintf(fmt{ic},real(val.(names{ic})(k,:)));
                rowStr=[rowStr str];
            end;
            mprintf(f,[rowStr '\n']);
            %             Write the complex parts
            if any(isComplex)
                cRowStr=[blanks(indent+iquant+indWidth) '+i'];
                for ic=1:nf2
                    nullStr=[blanks(widths(ic)-5) '--' blanks(3)];
                    if isComplex(ic)
                        str= sprintf(fmt{ic},imag(val.(names{ic})(k,:)));
                    else
                        str= repmat(nullStr,1,nCols(ic));
                    end;
                    cRowStr=[cRowStr str];
                end;
                mprintf(f,[cRowStr '\n']);
            end;
        end;
    end


    function WriteHex(val,rowSize)  % we assume val is an array of uint8
        indent=iquant;  % switch to a minimal indent
        ne=numel(val);
        nrows=ceil(ne/rowSize);
        nex=rowSize*nrows;
        val(nex)=0;  % extend the array
        %         val=reshape(val,hexRowSize,nrows)';
        pth=1;
        while pth<=ne
            pth1=min(ne,pth+rowSize-1);
            str=sprintf('%02X',val(pth:pth1));
            pth=pth1+1;
            %             str=sprintf('%02X',val(i,:));
            mprintf(f,'%s%s\n',blanks(indent),str);
        end;
        
    end

    function str=SizeString(sz,nf)
        %         Create a string like [1x2;5] (if nf is given)
        %         or [1x2] if nf is not given.  Here the size vector sz is assumed
        %         to be [1 2].
        str=sprintf('%dx',sz);
        str(end)=[];  % remove the x
        if nargin>1
            str=[str sprintf(';%d',nf)];
        end;
        str=['[' str ']'];
    end


end