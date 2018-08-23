function [vol norm]=rsFourierReconstruction(angles, fprojs, fnorms, symmetry, k)
% function [vol fnorm]=rsFourierReconstruction(angles, fprojs, fnorms, symmetry, k)
% Oversampled Fourier reconstruction for RSC data.  Usage:
% Given projection angles (in degrees), fourier-transformed class means, and fourier norms
% (both with zero frequency in the center), return the real-space volume.
% Optionally the back-projected normalization norm is returned, which can be used
% to show Fourier coverage.
%
% Example: We assume angles (nang x 3, in degrees), real-space projections
% (n x n x nangs) and the corresponding stack of ctfs,
% symmetry = 4;     % C4
% k=.1;             % Wiener constant
% nangs=size(angles,1);
% fprojs=single(zeros(n,n,nangs));
%   % Take the FT after shifting the particle to the origin.
% for i=1:nangs
%       Note the shifting of the projections to the origin before fft!
%     fprojs(:,:,i)=ctfs(:,:,i).*fftshift(fftn(ifftshift(projs(:,:,i))));
%     fnorms(:,:,i)=ctfs(:,:,i).^2;
% end;
% vol=rsFourierReconstruction(angles,fprojs,fnorms,symmetry);
% 
angs=rsDegToEuler(angles);  % convert to Euler angles in radians.
[n ny nangs]=size(fprojs);
norm=(zeros(n,n,n));
rvol=(zeros(n,n,n));
        rvol=OversamplingReconstr2(fprojs,angs,symmetry);
        norm=OversamplingReconstr2(fnorms,angs,symmetry);
fmsk=fuzzymask(n,3,n*0.4,n*0.05);

vol=fftshift(real(ifftn( fftn(ifftshift(rvol)).*ifftshift(fmsk)./(k^2+fftn(ifftshift(norm))) )));
% vol=(real(ifftn( fftn((rvol))./(k^2+fftn(ifftshift(norm))) )));

