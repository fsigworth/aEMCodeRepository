function ScalePDBCoords(inName,outName,sclFactor);
% Copies a pdb file, only scaling the X,Y,Z coordinates by sclFactor.

% examples of the lines that are read:
%REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000      
%ATOM      1  N   LYS A 331     -66.884  50.336  53.891  1.00195.42           N  

fin = fopen(inName,'rt');
if fin == -1   % couldn't be opened....
	error('Input file not found.');
end
fout=fopen(outName,'wt');
if fout<0
    error('Output file couldn''t be created.');
end;
fprintf(fout,'REMARK ---Copy of %s scaled by %8.3f\n',inName,sclFactor);
while (~feof(fin))
    s=fgetl(fin);
%      disp(s);

    s1=s;
    if numel(s)>53 && any(strcmp(s(1:4),{'ATOM' 'HETA'}))
        coords=zeros(1,3);
        for i=1:3
            ind=23+8*i;
            coords(i)=str2double(s(ind:ind+7));
        end;
%         coords=str2double(reshape(s(31:54),8,3)')'; doesn't work.
        newText=sprintf('%8.3f',sclFactor*coords);
        s1(31:54)=newText;
    end;
    fprintf(fout,'%s\n',s1);
%     disp(s1);
end;
fclose(fin);
fclose(fout);

 