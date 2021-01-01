function WriteStarFile(blockNames,blockData,fileName,headerText,activeFlags)
% function WriteStarFile(blockNames,blockData,fileName,headerText,activeFlags)
% Write a .star file with multiple data blocks. Each block uses a _loop
% construct, even if there is only one value per field. blockNames and
% blockData are cell arrays, one for each block. The optional headerText is
% a character row vector written at the beginning of the file. The
% optional activeFlags is a cell array of boolean vectors, one cell per
% block. Each boolean vector must have either 0 or nl elements, where nl is
% the number of lines in the corresponding block. A false active element
% says that line of output will not be written.
% Note that the output from ReadStarFile can be given immediately to WriteStar
% file to produce an equivalent star file, the only difference is that
% comment (and header text) information is lost.
% 
% Here's an example of arguments for writing a Relion 3.1 star file with
% optics and particles data blocks:
%   headerText='\n# version 30001\n';
%   blockNames={'optics'; 'particles'}; % Could also be given as 'data_optics' etc.
%   blockData={opt; pts};
%   activeFlags={[]; []};

nBlocks=numel(blockNames);
if nargin<5
    activeFlags=cell(nBlocks,1);
end;

fStar=fopen(fileName,'wt');

if nargin>3
    fprintf(fStar,headerText);
end;

for i=1:nBlocks
    WriteStarFileStruct(blockData{i},blockNames{i},fStar,activeFlags{i});
end;

fclose(fStar);
