function WriteMiText(mi,filename)
    % Write out an mi structure as a text file.  By default, filename is
    % constructed from the mi structure.  If filename is an empty string,
    % the output will be written to the command window.

if nargin<2
    filename=[mi.basePath mi.infoPath mi.baseName 'mi.txt'];
end;
if numel(filename)>0
    f=fopen(filename,'w');
else
    f=1;  % write to command window.
end;

names0=fieldnames(mi);
nFields=numel(names0);

% Order the fields in the mi structure, and label the types
% 0: string; 1: cell array of strings;  2: numeric; 3: ctf;
% 4: vesicle; 5: particle; 6: mask
orderedNames={
    'version' 0; 'identifier' 2; 'baseFilename' 0;'originalBasePath' 0;
    'basePath' 0; 'moviePath' 0;'imagePath' 0; 'procPath' 0; 'infoPath' 0;
    'stackPath' 0; 'tempPath' 0; 'movieFilename' 0;'gainRefName' 0;
    'imageFilenames' 1;'imageSize' 2;'pixA' 2;'kV' 2; 'cpe' 2; 'camera' 2;
    'tiltAngles' 2;'tiltAxis' 2;
    'frameDose' 2;'frameSets' 2;'frameShifts' 2;'doses' 2;'weights' 2;
    'ctf' 3;'mergeMatrix' 2;
    'vesicleModel' 2;'vesicle' 4;'particle' 5;'quality' 0; 'noiseModelPars' 2;
    'noiseModelCode' 1; 'boxSize' 2;
    'mask' 6;'notes' 1; 'log' 1};

nOrderedNames=size(orderedNames,1);
nameOrder=zeros(nOrderedNames,1);
for i=1:nOrderedNames
    q=strcmp(orderedNames{i,1},names0);
    ptr=find(q,1);
    if numel(ptr)<1
        ptr=0;
    end;
    nameOrder(i)=ptr;
end;

% Any fields that are missing from ordered Names are deleted
nameOrder(nameOrder==0)=[];

% Any leftover fields are put at the end
j=numel(nameOrder);
for i=1:nFields
    if ~any(nameOrder==i)
        j=j+1;
        nameOrder(j)=i;
    end;
end;

mi=orderfields(mi,nameOrder);


% names=fieldnames(mi);

% justNames=strjust(char(names),'right');  % right-justified
%
%
% % Write out the contents
% outName=[mi.basePath mi.infoPath mi.baseFilename 'mi.txt'];
% f=fopen(outName,'w');
% % if f<0
% %     error(['Couldn''t create the file ' outName]);
% % end;
% f=1
%
% for i=1:nFields
%     %%
%     lblanks=blanks(3);
%     mblanks=blanks(3+size(justNames,2));
%     name=names{i};
%     jname=[lblanks justNames(i,:)];
%     val=mi.(name);
%     switch name
%         case 'ctf'
%             sz=size(val);
%             mprintf(f,'%s [%1ux%1u] :\n',jname,sz(1),sz(2));
%             fldNames=fieldnames(val(1));
%             ctNames=(char(fldNames));
%             for j=1:size(ctNames,1)
%                 mprintf(f,'    %s%s = ',mblanks,ctNames(j,:));
%                 for k=1:numel(val)
%                     mprintf(f,'%10g ',val(k).(fldNames{j}));
%                 end;
%                 mprintf(f,'\n');
%             end;
%             mprintf(f,'\n');
%         case 'vesicle' % we'll show columns for each field
%             %             Convert the ok values to binary.
%             fields=fieldnames(val);
%             ncols=numel(fields);
%             nrows=numel(val.(fields{1}));
%             mprintf(f,'%s [%1ux%1u] : ',jname,nrows,1);
%             for j=1:ncols
%                 mprintf(f,'%10s',fields{j});
%             end;
%             mprintf(f,'\n');
%
%             if nrows>0
%                 mat=zeros(nrows,ncols);
%                 val1=val;
%                 val1.ok=val.ok(:,4);
%                 for j=3:-1:1
%                     val1.ok=val1.ok*2+val.ok(:,j);
%                 end;
%
%                 for j=1:ncols
%                     q=val1.(fields{j});
%                     mat(1:numel(q),j)=q;
%                 end;
%                 WriteMatrix(f,jname,mat,1);
%             end;
%             mprintf(f,'\n');
%         case 'particle'
%             fldNames=fieldnames(val);
%             for j=1:numel(fldNames)
%                 WriteMatrix(f,['   ' name '.' fldNames{j}],val.(fldNames{j}),0);
%             end;
%             mprintf(f,'\n');
%         otherwise  % general string or numeric
%             switch class(val);
%                 case 'char'
%                     mprintf(f,'%s = %s\n',jname,val);
%                 case 'cell'  % A cell array of strings or matrices
%                     if numel(val)==0
%                         mprintf(f,'%s = []\n',jname);
%                     else
%                         if isa(val{1},'char')  % strings
%                             mprintf(f,'%s = %s\n',jname,val{1});
%                             for j=2:numel(val)
%                                 mprintf(f,'%s   %s\n',mblanks,val{j});
%                             end;
%                             mprintf(f,'\n');
%                         else                   % matrices
%                             sz=size(val);
%                             mprintf(f,'%s {%1ux%1u} :\n',jname,sz(1),sz(2));
%                             for j=1:numel(val)
%                                 WriteMatrix(f,jname,val{j},1);
%                             end;
%                             mprintf(f,'\n');
%                         end;
%                     end;
%                 case {'single' 'double'}
%                     sz=size(val);
%                     if numel(sz)<3  % regular matrix
%                         if sz(1)<2  % single line
%                             if sz(2)<2 % scalar
%                                 if sz(2)==0
%                                     mprintf(f,'%s = []\n',jname);
%                                 else
%                                     mprintf(f,'%s = %g\n',jname,val);
%                                 end;
%                             else    % row vector
%                                 mprintf(f,'%s = [ ',jname);
%                                 mprintf(f,'%g  ',val);
%                                 mprintf(f,']\n');
%                             end;
%                         else
%                             WriteMatrix(f,jname,val,0);
%                             mprintf(f,'\n');
%                         end;
%                     else  % 3d matrix
%                         mprintf(f,'%s [%1ux%1ux%u] :\n',jname,sz(1),sz(2),sz(3));
%
%                         for j=1:size(val,3)
%                             WriteMatrix(f,jname,val(:,:,j),1);
%                             mprintf(f,'\n');
%                         end;
%                     end;
%
%             end;
%     end;
% end;
%

% function WriteMiText(mi)
hexRowSize=60;
maxStringsInLine=4;
maxNumbersInLine=10;
neForRowIndex=100;
iquant=4;

f=1;

fields=fieldnames(mi);
for fIndex=1:numel(fields)
    WriteField(mi.(fields{fIndex}),fields{fIndex},0);
end;



    function WriteField(val,fieldName,indent)
        c=class(val);
        switch c
            case 'struct'
                if strcmp(fieldName,'vesicle')
                    bigFieldNames={'x';'y'; 'shiftX'; 'shiftY';'ok';'r';'s'};
                    bigFieldWidths=[11 6; 11 6; 8 3; 8 3; 3 0; 13 5; 13 5];
                    WriteStructAsColumns(val,fieldName,bigFieldNames,bigFieldWidths,indent);
                else
                    WriteStruct(val,fieldName,indent);
                end;
            case {'single' 'double' 'int16' 'int32' 'uint8'}
                ne=numel(val);
                if ne<maxNumbersInLine
                    if ne<1
                        str=[blanks(indent) fieldName ' = []\n'];
                    else
                        str=[blanks(indent) fieldName ' = ' num2str(val(:)') '\n'];
                    end;
                    mprintf(f,str);
                    return
                else
                    sz=size(val);
                    str=[blanks(indent) fieldName ' ' SizeString(sz) ':\n'];
                    mprintf(f,str);
                    if strcmp(fieldName,'data')
                        WriteHex(val,indent+iquant);
                    else
                        %                         Handle matrices
                        
                        switch fieldName
                            case 'mergeMatrix'
                                rowSize=9;
                            case 'picks'
                                rowSize=sz(2);
                                val=val';
                            otherwise
                                rowSize=10;
                        end;
                        nrows=ceil(ne/rowSize);
                        showInd=ne>=neForRowIndex;
                        istr='';
                        ind=0;
                        ptr=1;
                        while ptr<=ne
                            ptr1=min(ne,ptr+rowSize-1);
                            str=sprintf('%12.6g ',val(ptr:ptr1));
                            ptr=ptr1+1;
                            if showInd
                                ind=ind+1;
                                istr=sprintf('%4d ',ind);
                            end;
                            mprintf(f,'%s%s%s\n',blanks(indent+iquant),istr,str);
                        end;
                    end;
                end;
            case 'char'
                mprintf(f,[blanks(indent) fieldName ' = ' val '\n']);
            case 'cell'
                ne=numel(val);
                sz=size(val);
%                 if ne<maxStringsInLine
%                     mprintf(f,[blanks(indent) fieldName ' = ' cellstr(val) '\n']);
%                 else
                    mprintf(f,[blanks(indent) fieldName ' ' SizeString(sz) ' cell:\n']);
                    for irow=1:ne
                        str=[blanks(indent+iquant) char(val{irow})];
                        mprintf(f,'%s\n',str);
                    end;
%                 end;
        end
    end

    function WriteStruct(val,name,indent)
        sz=size(val);
        ne=numel(val);
        str=[blanks(indent) name ' ' SizeString(sz) ' struct:'];
        mprintf(f,'%s\n',str);
        if ne>1  % array of stucts: call recursively
            for i=1:ne
                indName=sprintf('%s(%d)',name,i);
                WriteStruct(val(i),indName,indent+iquant);
            end;
        else
            names=fieldnames(val);
            for i=1:numel(names)
                WriteField(val.(names{i}),names{i},indent+iquant);
            end;
        end;
    end

    function WriteStructAsColumns(val,name,bigFieldNames,bigFieldWidths,indent)
%         special for mi.vesicle; also handles complex values.
        nrows=size(val.(bigFieldNames{1,1}),1);
        str=[blanks(indent) name ' ' SizeString(nrows) ' struct:'];
        mprintf(f,'%s\n',str);
        %       Find all the struct names
        allNames=fieldnames(val);
        littleFields={};
        j=0;
        for i=1:numel(allNames)
            isBig=strcmp(allNames{i},bigFieldNames);
            if any(isBig) && numel(val.(allNames{i}))==0
                bigFieldNames(isBig)=[];
                bigFieldWidths(isBig,:)=[];
                isBig=false;
            end;
            if ~any(isBig)
                j=j+1;
                littleFields{j,1}=allNames{i};
            end;
        end;
        
        ngps=numel(bigFieldNames);
        %         Create the header for the various fields
        ncols=zeros(ngps,1);
        for i=1:ngps
            ncols(i)=size(val.(bigFieldNames{i}),2);
        end;
        nameRow='';
        cmplx=false(ngps,1);
        for i=1:ngps
            gpWidth=bigFieldWidths(i,1)*ncols(i);
            
            str=sprintf(['%-' num2str(gpWidth) 's '],...
                [bigFieldNames{i} ' [' num2str(ncols(i)) ']']);
            nameRow=[nameRow str];
            cmplx(i)=any(imag(val.(bigFieldNames{i})(:))~=0);
        end;
        mprintf(f,[blanks(indent+iquant+5) nameRow '\n']);
        
        for j=1:nrows
            rowStr=[blanks(indent+iquant) sprintf('%4d ',j)];
            for i=1:ngps
                sz=bigFieldWidths(i,:);
                if sz(2)>0
                    fstr=['%' num2str(sz(1)-1) '.' num2str(sz(2)) 'g '];
                else
                    fstr=['%' num2str(sz(1)-1) 'd '];
                end;
                numStr= sprintf(fstr,real(val.(bigFieldNames{i})(j,:)));
                rowStr=[rowStr numStr];
            end;
            mprintf(f,[rowStr '\n']);
            %             Write the complex parts
            if any(cmplx)
                cRowStr=blanks(indent+iquant+5+1);  % one extra space for 'i'
                for i=1:ngps
                    sz=bigFieldWidths(i,:);
                    if cmplx(i)
                        fstr=['%' num2str(sz(1)-2) '.' num2str(sz(2)) 'gi '];
                        numStr= sprintf(fstr,imag(val.(bigFieldNames{i})(j,:)));
                    else
                        numStr=blanks(sz(1)*ncols(i));
                    end;
                    cRowStr=[cRowStr numStr];
                end;
                mprintf(f,[cRowStr '\n']);
            end;
        end;
        %        Finally, write any other fields
        for i=1:numel(littleFields)
            WriteField(val.(littleFields{i}),littleFields{i},indent+iquant);
        end;
        
    end
%         colWidths=cell2mat(bigFieldInfo{:,2});
%         ncols=cell2mat(bigFieldInfo(:,3});
%         strLength=strblanks(indent) name ' ' sizeString ' struct:'];
%         mprintf(f,'%s\n',str);
%
%         fnames=bigFieldInfo{:,1};
%         fcols=cell2mat(bigFieldInfo(:,2));
%         fstring=fieldInfo{:,3};
%         fstr1=
%         allNames=fieldnames(val);
%         littleFields={};
%         j=0;
%         for i=1:numel(allNames)
%             if ~strncmp(allNames{i},bigFieldInfo{:,1})
%                 j=j+1;
%                 littleFields{j,1}=allNames{i};
%             end;
%         end;
% %         Display the little fields
%
%         nfields=size(bigFieldInfo,1);
% %         Find the widest row
%         for i=1:nfields
%             nrows=max(nrows,size(val.(bigFieldInfo{i,1}),1));
%         end;
%
%         nrows=size(val.(bigFieldInfo{1,1}),1);
%         fnames=bigFieldInfo{:,1};
%         fsize=cell2mat(bigFieldInfo(:,2));
%         fstring=fieldInfo({:,3});
% %         Construct the header string
%
%
%         sz=size(val);
%         ne=numel(val);
%         str=[blanks(indent) name ' ' SizeString(sz) ' struct:'];
%         mprintf(f,'%s\n',str);
%         if ne>1  % array of stucts: call recursively
%             for i=1:ne
%                 indName=sprintf('%s(%d)',name,i);
%                 WriteStruct(val(i),indName,indent+iquant);
%             end;
%         else
%             names=fieldnames(val);
%             for i=1:numel(names)
%                 WriteField(val.(names{i}),names{i},indent+iquant);
%             end;
%         end;
%     end

    function WriteHex(val,indent)  % we assume val is an array of uint8
        indent=iquant;  % switch to a minimal indent
        ne=numel(val);
        nrows=ceil(ne/hexRowSize);
        nex=hexRowSize*nrows;
        val(nex)=0;  % extend the array
        %         val=reshape(val,hexRowSize,nrows)';
        ptr=1;
        while ptr<=ne
            ptr1=min(ne,ptr+hexRowSize-1);
            str=sprintf('%02X',val(ptr:ptr1));
            ptr=ptr1+1;
            %             str=sprintf('%02X',val(i,:));
            mprintf(f,'%s%s\n',blanks(indent),str);
        end;
        
    end

    function str=SizeString(sz)
        str=sprintf('%dx',sz);
        str(end)=[];
        str=['[' str ']'];
    end

%
%
%     function WriteMatrix(f,jname,val,noName)
%         len=numel(jname);
%         if noName
%             mprintf(f,'%s   [ ', blanks(len));
%         else
%             mprintf(f,'%s = [ ',jname);
%         end;
%         if size(val,1)>0
%             mprintf(f,'%10g ',val(1,:));
%             mprintf(f,'\n');
%             for i=2:size(val,1)
%                 mprintf(f,blanks(len+5));
%                 mprintf(f,'%10g ',val(i,:));
%                 if i==size(val,1)
%                     mprintf(f,']');
%                 end;
%                 mprintf(f,'\n');
%             end;
%         else
%             mprintf(f,']\n');
%         end;
%     end
%

end

