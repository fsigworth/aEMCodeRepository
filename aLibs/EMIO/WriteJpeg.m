function msc=WriteJpeg(m,filename,clipThreshold,doWrite)
% function msc=WriteJpeg(m,filename,clipThreshold,doWrite)
% Autoscale the image m and write it out as a .jpg or other graphics image, with the
% file type set by the extension of filename, .jpg, .png, or .tif
% clipThreshold is the fraction of gray
% values that are clipped at black and white ends of the range; default is
% 1e-3. If clipThreshold>1, it is the number of sds to include. The
% optional output msc is a clipped, 0..255 image ready for display using
% imaga(). If the optional argument doWrite=0, msc is returned but nothing 
% is written.

fileTypes={'.jpg','.tif','.png'};
if nargin<3
    clipThreshold=1e-3;
end;
if nargin<4
    doWrite=1;
end;
msf=imscale(m,256,clipThreshold);
msr=rot90(uint8(msf));

[~,~,ex]=fileparts(filename);
if ~any(strcmp(ex,fileTypes))
        disp(['WriteJpeg: Unrecognized extension in ' filename]);
return
end;

if doWrite
    imwrite(msr,filename);
end;

if nargout>0
    msc=max(0,min(msf,255)); % clip the output image.
end;
