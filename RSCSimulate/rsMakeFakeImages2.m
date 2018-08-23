function [imgs, y0s, partAngles]=rsMakeFakeImages2(map,padSize,vR,sm)
% function [imgs, y0s, trueAngles]=rsMakeFakeImages2(map,cropSize,vR,sm)
% Given the simulation parameters structure sm with fields
%   angles ni x 3 in degrees, [0 betas gammas]
%   bobs   ni x 1 in voxels
%   clicks ni x 2 in pixels
%   isos   ni x 1  (0 = right side out; 1 = inside out)
%   mbnOffset 1 x 1 in pixels; displacement of mbn from map center
% create a stack of ni noiseless images obtained as projections of map,
% padded to padSize. We assume that the particle's embedding location in
% the sphere is determined by sm.angles but its actual orientation
% (partAngles) differs due to rocking and also alpha errors resulting from
% click errors.  The projections are taken using partAngles, where
% partAlpha=residual rotation from click error, plus 180 degrees when
% inside-out, partBeta=angs(2)+rocks(2); and partGamma=angs(3).  The
% projections are then shifted by clicks (mirrored in x).
% y0s are the apparent yClick values as would be returned by the particle
% extractor.
degRads=180/pi;
floorCz=1;  % limit angles when particle center crosses the pole
n=size(map,1);
% Clicks are the displacement of the click to the true particle position.
% Thus -clicks is the vector from the image center to the particle, as the
% image is assumed to be centered on the click.
betas=sm.angles(:,2);
clicks=sm.clicks;  % click errors for each image.
% sm.isos are booleans, true= inside-out
effR=vR+sign(.5-sm.isos).*(sm.bobs-sm.mbnOffset);  % effective radius of particle center.
y0=effR.*sind(betas);  % true projected y-position of the particle
% True alphas will be nonzero rotating the particles to upright
% position due to the click error.
cz=hypot(clicks(1),clicks(2)+y0);  % always positive
% cz(abs(cz)<floorCz)=floorCz;
partAngles=sm.angles;
% positive clickX yields negative alpha.  Iso rotates everything by 180 too. 
partAngles(:,1)=sm.angles(:,1)+sm.rocks(:,1)+...
    degRads*real(asin(-clicks(:,1)./cz))+sm.isos*180;  % get rid of imag part when argument > 1
partAngles(:,2)=sm.angles(:,2)+sm.rocks(:,2);
% True alpha and beta will be influenced by rocking
% positive rockAlpha gives particle rotated right.

imgs0=Crop(rsMakeTemplates(partAngles,map),padSize,1); % pad the projections

% Shift left and down by click x and y.
imgs=real(ifft2(fft2(imgs0).*FourierShift(padSize,-clicks)));
% y0s=clicks(:,2)+y0;
y0s=cz;
