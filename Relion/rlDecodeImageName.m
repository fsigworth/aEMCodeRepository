function [slice,fileName]=rlDecodeImageName(str)
% Parse the Relion image name.
% For exampe, if str contains
% 000021@Micrograph.mrc
% we return the numeric value 21 and the string Micrograph.mrc
%
p=strfind(str,'@');
if numel(p)<1
    error(['Relion image string lacks ''@'' :' str]);
end;
p=p(1);
if p<2 || p>numel(str)-1
    error(['Relion image string is poorly formatted: ' str]);
end;
slice=str2double(str(1:p(1)-1));
fileName=str(p(1)+1:end);
