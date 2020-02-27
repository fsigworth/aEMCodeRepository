function [blockNames,blockData,ok]=ReadStarFile2(name)
% function [blockNames,blockData,ok]=ReadStarFile(name)
% -------This uses Matlab string arrays, hopefully to speed things up.-----
%----------
% Read a Relion .star file and put its contents into two cell arrays.
% blockNames is a cell array of the blockNames.  For example here are
% values from reading a post_run.star file.
% 
% >> blockNames{1} =
%   1�12 char array
% data_general
% >> blockData{1} =
%   struct with fields:
%                rlnFinalResolution: 7.0800
%       rlnBfactorUsedForSharpening: -215.7768
%         rlnFittedSlopeGuinierPlot: -53.9442
% >> blockNames{3} =
%   1�12 char array
% data_guinier
%  >> blockData{3} =
%   struct with fields:
%             rlnResolutionSquared: [51�1 double]
%         rlnLogAmplitudesOriginal: [51�1 double]
%                  rlnParticleName: [51�1 cell] % cell array of strings

% cd('/Users/fred/EMWork/relion13_tutorial/betagal/PrecalculatedResults/Refine3D')
% 
% name='post_run1.star';
tic

blockNames={};
blockData={};

commentMarkers={'#'};

ok=exist(name,'file');
if ~ok
    return
end;

lines=splitlines(string(fileread(name)));
% delete blank lines
% lines(lines=="")=[];

% delete from comment symbol to eol
hasComment=contains(lines,"#");
lines(hasComment)=extractBefore(lines(hasComment),'#');
lines=strip(lines); % remove white space before and after

nLines=numel(lines);

blockNames=lines;

% 
% fi=fopen(name);
% 
% nLines=0;
% C=cell(1,1);
% 
% % Load the whole file into the cell array C, handling comments
% disp('loading...');
% while ~feof(fi)
%     line=fgetl(fi);
%     p=strfind(line,commentMarkers);
%     hasComment = numel(p)>0;
%     if hasComment
%         line(p(1):end)=[];
%     end;
%     % Lines that only contain comments are treated as blank lines.
%     nLines=nLines+1;  % Count this line
%     if numel(line)<1
%         C{nLines}={{}}; % Count as blank
%     else
%         C{nLines}=textscan(line,'%s');
%     end;
%     if mod(nLines,100000)==1
%         disp(nLines)
%         disp(line);
%     end;
% end;
% fclose(fi);
% %% ------------------
% % C is a cell array {1,1} (a single '%s' is picked up) containing a cell
% % array {nc,1} where nc is the number of tokens in the line.

nBlocks=0;
P=1;  % line pointer

disp('scanning...');
while P<=nLines % loop through all the entries
    % skip blank lines
    while strlength(lines(P))<1
        P=P+1;
        if P>nLines
            return  % exit the function.
        end;
    end;
    % Get the block name
    % It starts with 'data_'
    if startsWith(lines(P),"data_")  % data block
        nBlocks=nBlocks+1;
        blockNames{nBlocks,1}=char(lines(P));
    else
        disp([num2str(nLines) ' ' lines(P)]);
        error(['''data_'' expected at line ' num2str(nLines) ' . exiting.'])
    end;
    P=P+1;
    if P>nLines  % reached the end of the file
        error(['End of file. Expected data at line ' num2str(nLines)]);
    end;
    
    while P<=nLines && strlength(lines(P))<1
        P=P+1;
    end;
    if P>nLines  % reached the end of the file
        error(['End of file. Expected data at line ' num2str(nLines)]);
    end;
 
    loopMode=lines(P)=="loop_";    
    if loopMode
        P=P+1;
                % skip blank lines
        while P<nLines && strlength(lines(P))<1
            P=P+1;
        end;
        if P>nLines
            error(['End of file.  Expected field name at line ' num2str(nLines)]);
        end;
    end;
    
    %  pick up fieldnames
    nFields=0;
    fieldNames=strings(0,1);
    fieldVals=strings(0,1);
    while P<=nLines && startsWith(lines(P),"_") % begins with underscore
        nFields=nFields+1;
        txt=extractAfter(lines(P),"_");
        if loopMode
            fieldNames(nFields,1)=txt;
        else
            strs=split(txt);
            fieldNames(nFields,1)=strs(1);
            if numel(strs)<1
                strs(2)="";
            end;
            fieldVals(nFields,1)=strs(2);
        end;
        P=P+1;
    end;
    
%     if P>nLines
%         return  % no data follows.
%     end;
    
    if loopMode  % Now the values follow immediately after the fieldnames
        disp('copying field values...');
% ---this is the slowest step in the whole function! ----
        nRows=0;
        fieldVals=strings(0,nFields);
        while P<=nLines && ~startsWith(lines(P),"data_")
            nRows=nRows+1;
            fieldVals(nRows,:)=split(lines(P));
            P=P+1;
        end;
        disp('Converting field values...');
    end;
toc
    return

%     Convert fieldVals to numeric when possible
    q=struct;
    for i=1:nFields
        fn=fieldNames{i};
        numericFVs=str2double(fieldVals(:,i));
        if all(~isnan(numericFVs))
            q.(fn)=numericFVs;
        else
            if nRows==1
                q.(fn)=fieldVals{1,i}; % field is a string
            else
                q.(fn)=fieldVals(:,i); % field is a cell array
            end;
        end;
    end;
   
    blockData{nBlocks,1}=q;
    
end;
disp('done.');
return;

