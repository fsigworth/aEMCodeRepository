function [d, vsh]=rsVesicleRadialDensity(m,mi,vesIndex)
thk=50;  % angstroms of membrane thickness
n=size(m);
ds=mi.imageSize(1)/n(1);
p=[mi.vesicle.x(vesIndex) mi.vesicle.y(vesIndex)]/ds+1;
sh=FourierShift(n,-p);
msh=fftshift(real(ifftn(fftn(m).*sh)));
r0=ceil(1.2*(mi.vesicle.r(vesIndex,1)+thk/mi.pixA)/ds);
vsh=Crop(msh,2*r0);
d=Radial(vsh);