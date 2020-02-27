% GratingFresnel2.m
% Use the Fresnel propagator to model electron waves below the sample
nx=512; % dimensions of outputs. Assumed to be even.
nz=2048;
psi=zeros(nx,nz);
origZ=nz/8;
origX=nx/2+1;
dx=1; % x per unit
dz=1; % z per unit

psi(origX,origZ)=1;

fpsi0=ones(nx,1);

for zi=1:nz-origZ % zi value
    zin=nz-origZ-zi+1; % z index in array
    fpsi(:,zin)=ifftn(
    












fx=20; % wavelength of grating, in pixels
% fy=20; % beam wavelength in pixels
amp=1; % extent of grating phase shift, radians





gratingY=.95; % height of grating in figure
gratingX=.2; % half width of grating in x

xs=(-nx/2:nx/2-1)'; % x coordinates; zero value is in the center
xss=fftshift(xs);
[mr, mi]=ndgrid(xs,ny:-1:0); % y coordinates increase downward
ny0=round(gratingY*ny)+1;  % y index of the grating

p=abs(xs)<gratingX*nx;  % define the width of the grating
phi1D=zeros(nx,1); % 1D phase perturbation
phi1D(p)=amp*(cos(xs(p)*2*pi/fx));
phi=zeros(nx,ny);
phi(:,ny0)=phi1D; % phase perturbation is only in one slice.

psi=zeros(nx,ny); % wave function
psi(:,ny)=1;  % top row is zero phase
psi0=psi;   % wave function with no perturbation

Psi=zeros(nx,ny); % FT of wave function
Psi(:,ny)=fft(psi(:,ny)); % top row
Psi0=zeros(nx,ny); % get the free propagation too.

dfz=4e-6; % scaling for diffraction angles

for i=ny-1:-1:1 % propagate waves downwardstep downward
    Psi(:,i)=fft(psi(:,i+1).*exp(1i*phi(:,i))) .* (exp(-1i* (pi*dfz*(xss.^2) + 2*pi/fy)));
    Psi0(:,i)=fft(psi0(:,i+1)).*fftshift(exp(-1i*(pi*dfz*(xs.^2) + 2*pi/fy)));
     psi(:,i)=ifft(Psi(:,i));
     psi0(:,i)=ifft(Psi0(:,i));
end;
%%
% psi(:,n0+1:n)=psi0(:,n0+1:n);
mysubplot(131);
imacsx(psi);
mysubplot(132);
imacsx(psi./psi0 - 1);
% imacsx(psi./psi0);
mysubplot(133);
imags(abs(psi))

