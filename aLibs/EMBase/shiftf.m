function out=shiftf(in,delta)
% function out=shiftf(in,delta)
% Analogous to circshift, but uses the fft to  
% perform fractional shifts of the image or image stack in.
% delta = [dx dy] is the shift, with multiple rows for multiple images.
[nx, ny, nim]=size(in);
[X,Y]=ndgrid((-nx/2:nx/2-1)/nx,(-ny/2:ny/2-1)/ny);
out=zeros([nx ny nim],'single');
for i=1:nim
    P=exp(1j*2*pi*fftshift((delta(i,1)*X+delta(i,2)*Y)));
    out(:,:,i)=real(ifftn(fftn(in(:,:,i)).*P));
end
