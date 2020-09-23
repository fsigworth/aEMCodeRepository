function msc=WriteJpeg(m,filename,clipThreshold,doWrite)
% function msc=WriteJpeg(m,filename,clipThreshold,doWrite)
% Autoscale the image m and write it out as a jpeg image.  If filename has
% no extension add '.jpg' to it.  clipThreshold is the fraction of gray
% values that are clipped at black and white ends of the range; default is
% 1e-3. If clipThreshold>1, it is the number of sds to include. The
% optional output msc is a clipped, 0..255 image ready for display using
% imaga(). If the optional argument doWrite=0, msc is returned but nothing 
% is written.

if nargin<3
    clipThreshold=1e-3;
end;
if nargin<4
    doWrite=1;
end;
if numel(regexp(filename,'.+\.jpg'))+numel(regexp(filename,'.+\.jpeg'))==0
    filename=[filename '.jpg']; % add .jpg if not present.
end;
msc=imscale(m,256,clipThreshold);
if doWrite
    imwrite(rot90(uint8(msc)),filename);
end;

if nargout>0
    msc=max(0,min(msc,255)); % clip the output image.
end;
