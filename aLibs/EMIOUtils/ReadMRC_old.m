function [map,s]=ReadMRC(filename);
% function [map,s]=ReadMRC(filename);
% Read a 3D map from an .mrc or CCP4 map file
% and return a structure containing various parameters read from the file.
% This function reads 2d and 3d real maps in byte, int16 and float32 modes.
% Modified 25 Jan 06 to read both little and big endian files.

% We initially specify little-endian data:
f=fopen(filename,'r','ieee-be');

% Get the first 10 values, which are integers:
% nc nr ns mode ncstart nrstart nsstart nx ny nz
[a,cnt]=fread(f,10,'int32');

% a(1) is the number of columns.  See if the number is reasonable.
if abs(a(1))>1e5  % must be big-endian! Go open the file again.
  fclose(f);
  f=fopen(filename,'r','ieee-le');
  [a,cnt]=fread(f,10,'int32');
end;

% a(1:10)

mode=a(4);
% mode

% Get the next 12, which are floats.
% the first three are the cell dimensions.
% xlength ylength zlength alpha beta gamma mapc mapr maps
% amin amax amean.
[b,cnt]=fread(f,12,'float32');
%  b
% rez=b(1)/a(1);
s.rez=b(1)/a(1);
s.celldim=b(1:3);

% get the next 30, which brings us up to entry no. 52.
[c,cnt]=fread(f,30,'int32');
% c(1:3)

% the next two are supposed to be character strings.
[d,cnt]=fread(f,8,'char');
% d
% char(d(1:3))
% char(d(5:8))

% Two more ints...
[e,cnt]=fread(f,2,'int32');

% up to 10 strings....
ns=min(e(2),10);
for i=1:10
  [g,cnt]=fread(f,80,'char');
  str(i,:)=char(g)';
end;

% disp('header:'); disp(' ');
% disp(str(1:ns,:));
% disp(' ');

s.header=str(1:ns,:);

% Get ready to read the data.
s.nx=a(1);
s.ny=a(2);
s.nz=a(3);

switch mode
  case 0
    string='int8';
  case 1
    string='int16';
  case 2
    string='float32';
  otherwise
    disp('ReadMRC: unknown data mode');
    mode
    string='???';
end;

[map,cnt]=fread(f,s.nx*s.ny*s.nz,string);

map=reshape(map,s.nx,s.ny,s.nz);
fclose(f);
