function [slice,fileName]=rlDecodeImageName(str)
% function [slice,fileName]=rlDecodeImageName(str)
% function [sliceNums,fileNames]=rlDecodeImageName(imgNames)
%
% Parse the Relion image name as found in a Relion star file.
% For exampe, if str contains
% 000021@Micrograph.mrc
% we return the numeric value 21 and the string Micrograph.mrc
% If the input is a cell array of strings, returned is a vector sliceNums
% and a cell array of strings fileNames.
%
% Example: put the image number and stack name into new fields of the
% struct.
% [names,data]=ReadStarFile('file.star');
% q=data{1};
% [q.rlnImageNos,q.rlnStackNames]=rlDecodeImageName(q.rlnImageName);

if isa(str,'cell')
    inp=str;
    ns=numel(inp);
    slice=zeros(ns,1);
    fileName=cell(ns,1);
else
    inp={str};
    ns=1;
    slice=0;
end;
for i=1:ns
    str=inp{i};
    p=strfind(str,'@');
    if numel(p)<1
        error(['Relion image string lacks ''@'' :' str]);
    end;
    p=p(1);
    if p<2 || p>numel(str)-1
        error(['Relion image string is poorly formatted: ' str]);
    end;
    slice(i)=str2double(str(1:p(1)-1));
    name=str(p(1)+1:end);
    if ns>1
        fileName{i,1}=name;
    else
        fileName=name;
    end;
end;
