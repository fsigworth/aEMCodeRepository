function names=deGetRawImageNames(folderName)
% Create a cell array of filenames for DE-12 raw images, given the folder that contains them.
% We assume that filenames are of the form RawImage_N.tif where N is
% numeric.
d=dir(folderName);
if numel(d)<3
    error(['No files found in this path: ' folderName]);
end;
names={};

% Scan and read image numbers
for i=3:numel(d)
    name=d(i).name;
    q=strfind(name,'.tif');  % it's a tif file
    if numel(q)>0 && q(1)>1
        nm=name(1:q(1)-1);
        q1=strfind(name,'RawImage_');
        if numel(q1)>0       % it's a raw image
            imgNumber=str2double(nm(q1(1)+9:numel(nm)))+1;
            names{imgNumber}=[folderName d(i).name];
        end;
    end;
end;
% % If any names are missing or blank, delete those entries from the array.
j=1;
while j<=numel(names)
    if numel(names{j})<1 % no name there
        names(j)=[];
    else
        j=j+1;
    end;
end;
