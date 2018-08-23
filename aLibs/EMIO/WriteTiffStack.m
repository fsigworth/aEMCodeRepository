function WriteTiffStack(m,pixA,filename)
% function WriteTiffStack(m,pixA,filename)
% Write out the uint8 stack of images m to an LZW-compressed TIFF file.
% If an extension is lacking, append '.tif'
if numel(regexp(filename,'.+\.tif'))+numel(regexp(filename,'.+\.tiff'))==0
    filename=[filename '.tif'];
end;
if ~isa(m,'uint8')
    m=uint8(m);
end;
[nx,ny,nz]=size(m);

% We'll record the pixel size in the ImageDescription string.
str=sprintf('s.pixA= %f; s.nx= %d; s.ny=%d; s.nz=%d;',pixA,nx,ny,nz);
% Write the first slice
imwrite(m(:,:,1),filename,'compression','lzw','description',str);
for i=2:nz  % append the rest
    imwrite(m(:,:,i),filename,'compression','lzw','writemode','append');
end;

