function [grid,hole,shot]=EPUDecodeName(name)
% Given an EPU name such as
%   Info_m/GridSquare_28359688_FoilHole_28373575_01mi.txt
% return the numberic values 28359688, 28373575, 1

% Skip any path specifications
p=strfind(name,filesep);
if numel(p)>0
    name=name(p(end)+1:end);
end;

% pick up alternating strings and numbers
C=textscan(name,'%s %d %s %d %d','Delimiter','_');
% celldisp(C)
grid=C{2};
hole=C{4};
shot=C{5};
