function [m,s,ok,rawImgs]=ReadZTiff(filename,startSlice,numSlices)
% Read a compressed TIFF file using the NNR compression.  The output m is
% in singles.  It stops with an error if the file doesn't contain the
% compression data; alternatively, if the ok output argument is given, it
% returns ok=false in that case without giving an error. The returned
% struct s contains the pixA field, i.e. angstroms per pixel, along with a
% size field [nx ny nframes], also redundant fields nx, ny, nz, and
% various other information including the compression parameters.
% If numSlices=0 then m is returned as [] but s has the information.

% Here is an example of the contents of the idString that is read from the file.  
% (There are no leading spaces or % in the lines of the string, those are
% here only for formatting):
%   s.compVersion=1.01;
%   s.origPixA=1.3489;
%   s.size=[4096  4096     1];
%   s.origSize=[4096  4096];
%   s.origClass='single';
%   s.fs=[0.07            0.1           0.14            0.2            0.4];
%   s.as=[22534.014     -1786.4751      4701.7706      327.79427      156.54706];
% s.size is the size of the file to be returned [nx ny nframes].
% s.origSize is the same as s.size unless we have downsampled the original
% image to save more space.  (This feature is not yet implemented.)
% s.origClass is the data type to be returned, e.g. 'single', 'int16', 'uint8'
% s.fs are the standard-deviations of Gaussians in frequency. In the scaling
% used, 0.5 is Nyquist.
% s.as are amplitudes of the Gaussians obtained in the fit to the original
% power spectrum.  These are scaled such that a filter using sqrt(s.as) as
% coefficients will restore the original image's amplitude.

if nargin<3
    numSlices=inf;
end;
if nargin<2
    startSlice=1;
end;
startSlice=max(1,startSlice);

saveRawImgs=nargout>3;

ok=true;
m=[];
s=struct;

try
t=Tiff(filename,'r');  % open the file and get the TIFF object t.
catch
    if nargout>2
        ok=false;
        return
    else
        error(['Could not open TIFF file: ' filename]);
    end;
end;

% Count the number of pages in the file
lastSlice=1;
while ~t.lastDirectory
    t.nextDirectory;
    lastSlice=lastSlice+1;
end;
s.nz=lastSlice;  % in case we won't get it from reading

% Cycle through the file pages
if startSlice>lastSlice
    error('startSlice is larger than file length');
end;
numSlices=min(numSlices,lastSlice-startSlice+1);
endSlice=startSlice+max(0,numSlices-1);  % read from at least one.

for iimg=startSlice:endSlice
    t.setDirectory(double(iimg)); % requires double argument!
    %     Read the parameters
    try
        idString=t.getTag('ImageDescription');
    catch
        if nargout>2
            ok=false;
            return
        else
            error('Not a z-compressed tiff file');
        end;
    end;
    eval(idString);  % fill in fields of s
    if ~isfield(s,'compVersion')
        if nargout>2
            ok=false;
            return
        else
            error('Not a z-compressed tiff file');
        end;
    end;

    %   -------Read and image and do the decompression------
    if numSlices>0  % if we do have something to read
        crf=t.read;
        if iimg==startSlice;  % first one: allocate output of the desired
%                               numeric type.
            m=zeros([s.origSize(1:2) endSlice-startSlice+1],s.origClass);
            if saveRawImgs
                rawImgs=zeros(size(m),'uint8');
            end;
        end;
        if saveRawImgs
            rawImgs(:,:,iimg-startSlice+1)=crf;
        end;
        rf=double(crf)-128;  % reconstruct using double-precision initially
        s.nx=size(rf,1);
        s.ny=size(rf,2);
        fsModel=ccdAliasedGaussians(s.nx,s.fs,s.as);
        m(:,:,iimg-startSlice+1)=real(ifftn(fftn(rf).*ifftshift(sqrt(fsModel))));
    elseif s.compVersion<1.01 % We have to assume size = origSize
        s.nx=s.origSize(1);
        s.ny=s.origSize(2);
    else  %  We can read the size variable.
        s.nx=s.size(1);
        s.ny=s.size(2);
    end; % if nSlices>0
end; 
close(t);
s.pixA=s.origPixA*s.origSize(1)/s.nx;

