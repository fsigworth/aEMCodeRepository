
function [m, pixA, ok]=ReadEMFile(name,startSlice,numSlices)
% function [m, pixA, ok]=ReadEMFile(name)
% function m=ReadEMFile;   % simplest form: puts up a file selector
% --General image file reader--
% Based on the file extension, read the image and return the pixel size
% (in angstroms).  The returned value of pixA is zero if a pixel size is
% not available.
% Note that the original data type (uint8, uint16, single) is returned!  No
% conversions are done.  However, rotations are made to standard graphics
% files to correspond to the Cartesian coordinate convention.
% This function understands the following extensions:
% z.tif (compressed Tiff)
% .tif  (generic, or TVIPS files) --rotated 270 degrees
% .dm3, .dm4
% .mrc
% .hed, .img
% .jpg   --rotated 270 degrees
% The entire image is returned for .dm3, .dm4 and .jpg files, regardless of
% the startSlice and numSlices parameters. for TVIPS .tif files, either all
% data (numSlices>0) or no data, just s, are returned.
% If the file has the wrong extension, m is returned as [] and ok=0.
% fs Jan 2011, Aug 2013
% Added check for existence May 2014.
% Added startSlice, numSlices, ZTiff support Feb 2015
% 

pixA=0;
m=[];
ok=0;
showInfo=false;

if nargin<1
    fprintf('Opening the file: ');
    typeStr='*.mrc;*.mrcs;*.st;*.dm3;*.dm4;*.tif;*.hed;*.img;*.jpg';
    [name,path]=uigetfile(typeStr,'Open a file');
    if isnumeric(name)
        disp(' no file selected.');
        return
    else
        disp([path name]);
        showInfo=true;
    end;
    cd(path);
end;
if nargin<2
    startSlice=1;
end;
if nargin<3
    numSlices=inf;
end;
someData=numSlices>0;



name=InterpretFilename(name);
if ~exist(name,'file')
    return
end;
[pa, nm, ex]=fileparts(name);
isZtif=strcmp(name(end-4:end),'z.tif');
pixA=0;  % default value
m=zeros(100,100); % default value
ok=1;
switch lower(ex)
    case {'.mrc','.mrcs','.st'}
        [m, s]=ReadMRC(name,startSlice,numSlices);
        pixA=s.rez/s.nx;
    case {'.dm3','.dm4'}
        [m, pixnm]=ReadDMFile(name);  % Can't read just part of the file
        pixA=10*pixnm;
    case '.tif'
        if isZtif
            [m,s,ok]=ReadZTiff(name,startSlice,numSlices);
            if ok 
                pixA=s.pixA;
            end;
        end;
        if ~isZtif || ~ok  % try for a TVIPS Tiff file
            s=ReadTiffMetadata(name);  % Get metadata from a TVIPS file.
            pixA=s.pixA;               % will be zero if not a TVIPs file
            %         Although TIFF doesn't support it, TVIPS 2-byte files are encoded
            %         as signed integers, so we coerce the type.
            if someData
                m=rot90(imread(name),3);
                if isa(m,'uint16')&&(pixA>0)
                    sz=size(m);
                    m=typecast(m(:),'int16');
                    m=reshape(m,sz);
                end;
            end;
        end;
    case  {'.hed','.img'}
        m=ReadImagic(name,startSlice,numSlices);
    case '.jpg'
        m=rot90(imread(name),3);  % allow jpegs to be read too.
    otherwise
        warning(['Unknown file type: ' ex]);
        ok=0;
end;
if showInfo && numel(m)>0
    whos m
end;