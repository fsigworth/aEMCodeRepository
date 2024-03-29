function [yest, indices, aldata] = gridAlign( data, reference, nrot, symmetry )
% arAlign:  re-estimate the model square matrix yest using the maxima
% of the cross-correlation with the reference matrix.

[n n1 m] = size(data);
if nargout==3
  aldata=zeros(n,n,m);
end;

rotsums = zeros(n,n,nrot);
qrot = 2*pi / (nrot*symmetry);	% quantum of rotation

% Make rotated copies of the reference
for r=1:nrot
  rpre = grotate( reference, (r-1)*qrot );
  fpre(:,:,r) = conj(fft2( rpre )); % We will use this to do correlation.
imacs(rpre); drawnow; title(r);
end

for i=1:m

  % Do the correlation
  fdat=fft2(data(:,:,i)); % FT of datum
  for r = 1:nrot
    ycorr(:,:,r) = real(ifft2(fdat.*fpre(:,:,r) )); % cross-correlation of ypre and datum
  end
  [maxval, is, js, rs] = max3d(ycorr);

  % store the alignment parameters
  indices(1:3,i)=[is js rs]';

  if nargout<3
    % we only want the mean; do the translational alignment
    rotsums(:,:,rs) = rotsums(:,:,rs) + circshift( data(:,:,i), 1-is, 1-js );
  else
    aldata(:,:,i)=grotate(data(:,:,i),(1-r)*qrot);
imacs(aldata(:,:,i)); drawnow; title(i);
  end;
end

if nargout<3
  yest = zeros(n);
  for r = 1:nrot
    yest=yest + mrotate( fest(:,:,r), (1-r)*qrot );	% Rotate each reconstruction back.
  end
else
  yest=squeeze(sum(aldata,3));
  yest=yest/m;
end;
