function [imgs y0s]=rsMakeFakeImages(map,cropSize,vR,angs,rocks,bobs,clicks)
% Given the ni x 3 angles [0 betas gammas] and ni x 2 rock angles 
% [alphas betas], along with vectors of bobs and ni x 2 clicks (both in
% pixels) create a stack of
% ni noiseless images obtained as projections of map.
% Images are rotated according to
%   alpha=residual rotation from click error
%   beta=angs(2)+rocks(2);
%   gamma=angs(3)
% and then shifted by clicks (mirrored in x).

degRads=180/pi;
floorCz=0.1;
n=size(map,1);
if cropSize==0
    cropSize=n;
end;
% clicks are the displacement of the click to the true particle position.
% Thus -clicks is the vector from the image center to the particle after
% centering the image on the click.
betas=angs(:,2);
y0=(vR+bobs).*sind(betas);  % true projected y-position of the particle
    % Calculate the alpha residual due to the click error
cz=hypot(-clicks(1),-clicks(2)+y0);
cz(abs(cz)<floorCz)=floorCz;
angs(:,1)=degRads*real(asin(-clicks(:,1)./cz));  % get rid of imag part when argument > 1
angs(:,1:2)=angs(:,1:2)+rocks;

imgs0=Crop(rsMakeTemplates(angs,map),cropSize,1);

imgs=real(ifft2(fft2(imgs0).*FourierShift(cropSize,-clicks)));
y0s=clicks(:,2)+y0;