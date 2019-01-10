% ConvertMRCToTiff.m
% Convert an .mrc file to .tif
% read an mrc file in any format, and write a corresponding 8 or 16-bit
% TIFF file.
% We handles modes 0 and 1, which are signed integers, by converting them
% to unsigned for the TIFF format.  We handle mode 2, which is floating
% point, by scaling it up if its span is less than the span variable.
% Modes 6 and 101 are assumed to give unsigned integer data, ready for TIFF.

span=1000;  % minimum span for the integer conversion of a floating-point file.

% if desired, put up a file selector...
[name,path]=uigetfile('*.mrc');
cd(path);
%  ...alternatively, just assign the name variable:
% name='myfile.mrc';

% Read the mrc file
[m,s]=ReadMRC(name);

m=single(m); % convert to floating point
mi=min(mf(:));
mx=max(mf(:));
if mi<0
    m=m-mi; % force it to be non-negative
    mx=mx-mi; % correct the maximum
end;

switch s.mode
    case {0 101}
        m=uint8(m); % force it to be unsigned.
    case {1 6} % 16 bit
        m=uint16(m); % force it to be unsigned, as TIFF wants.
    case 2  % it's floating point, expand the range to span
        if mx<span
            m=m*span/abs(mx);
        end;
end;

% create a new filename that ends in .tif
[~,baseName,ext]=fileparts(name);
newName=[baseName '.tif'];

% Write the TIFF file
imwrite(m,newName,'tiff');

disp([newName ' written.']);

