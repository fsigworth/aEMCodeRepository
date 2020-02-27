function handle=WriteMRCHeader(map,pixA,filename,sizes,org,mode,mz)
% function handle=WriteMRCHeader(map,pixA,filename,sizes,org,mode,mz)
% Write out the header of an MRC map file, and leave the file open,
% returning the file handle, so that data can be written sequentially into
% the file and then the file closed.  Data are written in little-endian style,
% with the datatype specified by mode.  By default mode=2, which means
% float32 values.  Other options are 0: uint8; 1: int16; 32: packed uint4.
% If you want to write more slices than are contained in map, give the
% total number of slices (or 2d images) in the optional sizes argument.
% Alternatively you can give the total size of the map to be written as
% sizes=[nx ny nz]. The optional argument org is the starting point of the
% volume, an integer in units of voxels.  The default value is [c c c]
% where c=-floor(n/2). The argument mz is the number of planes per volume.
% Thus the number of volumes in the file is sizes(3)/mz.  For images, mz=1
% (default).
%
% Example: write out 10,000 random images.
%     images=randn(64,64,1000);
%     f=WriteMRCHeader(images,2.8,'test.mrc',10000);
%     fwrite(f,images,'float32');
%     for i=2:10
%         images=randn(64,64,1000);
%         fwrite(f,images,'float32');
%     end;
%     fclose(f);
% 
% Changed Jul2011 to support the origin argument.  fs
% Changed Dec2011 to support the mode argument.
% Changed Sep2013 to support 3-element sizes argument.
% Files are always written in little-ended format.
% Figure out if we have a little-ended machine.

q=typecast(int32(1),'uint8');
machineLE=(q(1)==1);  % true for little-endian machine

hdr=int32(zeros(256,1));
if nargin<3
    sizes=size(map);
elseif numel(sizes)==1  % only gives the number of planes
    sizes=[size(map,1) size(map,2) sizes];
elseif numel(sizes)==2 % only one plane
    sizes(3)=1;
end;

if nargin<5
    org=-floor(sizes/2);  % Default origin is such that 
else
    org=org(:);
end;
if nargin<6
    mode=2;
end;
if nargin<7
    mz=1;  % by default, we are a stack of images.
end;
%     
% Get statistics.
if mode==32
    corr=1/16;  % approximate correction for packed data
else
    corr=1;
end;
map=reshape(map,numel(map),1);  % convert it into a 1D vector
% theMean=mean(single(map))*corr;
% theSD=std(single(map))*corr;
% theMax=single(max(map))*corr;
% theMin=single(min(map))*corr;
theMean=0;
theSD=0;
theMax=0;
theMin=0;

% Number of voxels in the map or submap
hdr(1:3)=sizes; % number of columns, rows, sections
hdr(4)=mode;  % mode: real, float values
hdr(5:7)=org;  % origin of the map

% Number of voxels in the full unit cell
hdr(8:10)=[hdr(1:2)' mz];  % number of intervals along x,y,z
% Unit cell sizes and angles
hdr(11:13)=typecast(single(single(hdr(8:10))*pixA),'int32');  % Cell dimensions
hdr(14:16)=typecast(single([90 90 90]),'int32');   % Angles
hdr(17:19)=(1:3)';  % Axis assignments

hdr(20:22)=typecast(single([theMin theMax theMean]'),'int32');
hdr(23)=0;  % Space group 0 (default)
if machineLE
    hdr(53)=typecast(uint8('MAP '),'int32');
    hdr(54)=typecast(uint8([68 65 0 0]),'int32');  % LE machine stamp.
else
    hdr(53)=typecast(uint8(' PAM'),'int32');  % LE machine stamp, for writing with BE machine.
    hdr(54)=typecast(uint8([0 0 65 68]),'int32');
end

hdr(55)=typecast(single(theSD),'int32');

handle=fopen(filename,'w','ieee-le');  % Here we force little-ended order.
if handle<0
    error(['File could not be opened for writing: ' filename]);
end;

count1=fwrite(handle,hdr,'int32');
