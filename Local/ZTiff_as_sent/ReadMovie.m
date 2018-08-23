function [m,s,ok]=ReadMovie(name,start,num)
% function [m,s]=ReadMovie(name,start,num)
% File reader for k2 or f2 camera movies, .mrc, .tif, z.tif.
% Based on the file extension, read the image and return the struct s
% which has fields nx,ny,nz and, in some cases, pixA.  When no pixA value
% is available, as in generic TIFF files, s.pixA is returned as zero.
% Also handles gzipped files with .gz following the true extension.
if nargin<2
    start=1;
end;
if nargin<3
    num=inf;
end;
m=0; % default values
s=struct;
ok=true;
[~,~,ex]=fileparts(name);
% % Special case for gzipped files:
% isGzFile=strcmp(ex,'.gz');
% if isGzFile
%     gzName=name;
%     name=[AddSlash(pa) nm];
%     [~,~,ex]=fileparts(name);
%     if ~exist(name,'file')  % not already unzipped
%         zipErr=system(['gunzip -kf ' gzName]);
%         if zipErr
%             ok=false;
%         return
%         end;
%     end;
% end;

% Special case of ztiff files
if strcmp(name(end-4:end),'z.tif')
    ex='z.tif';
end;

switch lower(ex)
    case {'.mrc','.st','.mrcs'}
        [m, s]=ReadMRC(name,start,num);
    case '.tif'
        [m, s]=ReadTiffStack(name,start,num);
    case 'z.tif'
        [m,s,ok]=ReadZTiff(name,start,num);
    otherwise
        if nargout<3
            error(['Unknown file type: ' ex]);
        else
            ok=false;
        end;
end;
% if isGzFile % We'll delete the unzipped file
%     disp(['deleting ' name]);
% %     delete(name);
% end;
