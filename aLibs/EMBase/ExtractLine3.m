function line = ExtractLine3(m, n, vector, origin)
% function line = ExtractLine3(m, n, vector, origin)
% Extract, as an n x 1 vector, values of the 3D volume m along a line
% of length n, with the center of the line from map(origin) and with the
% points of the line beyond the center being in the direction of vector
% from the origin. (vector and origin are 1 x 3 vectors.) Trilinear interpolation is used. Based on the code of
% ERotate3. If coordinates go out of bounds, we clip to the nearest
% in-bound coordinate.

nearest=0; % nearest-neighbor interpolation instead of trilinear.

[nx,ny,nz]=size(m);

if nargin<4 || numel(origin)==0
    origin=ceil((size(m)+1)/2);  % fft center
end;

% normalize the direction vector
vNorm=vector(:)'/sqrt(vector(:)'*vector(:));
vSteps=ceil(-n/2):floor((n-1)/2);
rampV=origin+vSteps'*vNorm;

% Convert to 1D and clip the coordinates to be in bounds.
epsl=1e-6; % We insure that rounding up and down leaves us in bounds.
epsh=1+1e-6;
X=max(min(rampV(:,1),nx-epsh),1+epsl);
Y=max(min(rampV(:,2),ny-epsh),1+epsl);
Z=max(min(rampV(:,3),nz-epsh),1+epsl);

if nearest
    % No interpolation:
    line= m(round(X)+nx*(round(Y-1)+ny*round(Z-1)));
else
    % Linear interpolation:
    X0=floor(X);
    Xi=X-X0;
    X1=X0+1;
    
    Y0=floor(Y);
    Yi=Y-Y0;
    Y1=Y0+1;
    Y0=(Y0-1)*nx;  % Convert to linear index component for rows.
    Y1=(Y1-1)*nx;
    
    Z0=floor(Z);
    Zi=Z-Z0;
    Z1=Z0+1;
    Z0=(Z0-1)*nx*ny;  % Convert to linear index component for planes.
    Z1=(Z1-1)*nx*ny;
    
    line=(1-Zi).*((1-Yi).*((1-Xi).*m(X0+Y0+Z0) + Xi.*m(X1+Y0+Z0))...
        + Yi.*((1-Xi).*m(X0+Y1+Z0) + Xi.*m(X1+Y1+Z0)))...
        +Zi.*((1-Yi).*((1-Xi).*m(X0+Y0+Z1) + Xi.*m(X1+Y0+Z1))...
        + Yi.*((1-Xi).*m(X0+Y1+Z1) + Xi.*m(X1+Y1+Z1)));
end;

