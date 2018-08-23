function [coords types]=ReadPDBAtoms2(filename, hetatms, chains)
% function [coords types]=ReadPDBAtoms(filename, hetatms);
% Fast pdb coordinate reader.  Reads the coordinates of all ATOM records or
% both ATOM and HETATM records (if hetatms=1, default).  Given na atoms,
% coords is a 3 x na array, and types is an na x 4 character array, with
% each type padded with trailing blanks.
% Code is modified from pdbread() in the Bioinformatics Toolbox.
% fs 26 Apr 11

if nargin<2
    hetatms=1;
end;

if nargin<3
    chains='';
    allChains=true;
else
    allChains=false;
end;

if hetatms>0
    atomText={'ATOM','HETATM'};
else
    atomText='ATOM';
end;

fid = fopen(filename,'r');

if fid == -1,
    error('ReadPDBAtoms:CouldNotOpenFile',...
        'Could not open file %s.', filename);
end;

theLines = textscan(fid,'%s','delimiter','\n');
theLines = theLines{1};  % textscan returns a 1x1 cell containing a cell array.
fclose(fid);

nl=numel(theLines);

% Count the number of atom coordinates
na=0;
okLines=false(1,nl);
for i=1:nl
    tline = theLines{i};
    if ischar(tline)
        len = length(tline);
        if len>54 && any(strcmpi(deblank(tline(1:6)),atomText)) && ...
                (allChains || numel(strfind(chains,tline(22)))>0)
            okLines(i)=true;
        end;
    end;
end;

coords=zeros(3,sum(okLines));
types=char(zeros(na,4));
na=0;
inds=find(okLines);
for i=inds
    tline = theLines{i};
%     len = length(tline0)
    
    %     % RCSB web site format requires each line to have 80 characters. This
    %     % avoids exceeding the matrix dimension for lines with
    %     % Less than 80 characters.
%     tline = [tline0 blanks(80-len)];
%     Record_name = deblank(tline(1:6));
    
%     switch upper(Record_name)
%         %
%         case atomText
%             if allChains || numel(strfind(chains,chain))>0
            na=na+1;
            s=strtrim(tline(13:16));  % atom type
            sl=length(s);
            types(na,:)=[s blanks(4-sl)];
            coords(1,na) = sscanf(tline(31:38),'%f');
            coords(2,na) = sscanf(tline(39:46),'%f');
            coords(3,na) = sscanf(tline(47:54),'%f');
%             end;
%     end;
end;

% plot3(coords(1,:),coords(2,:),coords(3,:),'k.');