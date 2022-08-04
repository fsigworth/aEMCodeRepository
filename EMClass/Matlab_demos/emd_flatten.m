% emd_flatten.m
% Manual solvent-flattening of the TRPV1 map

[m,s]=ReadMRC('emd_5778.map');
f0=.01;
H=1+1./(f.^2/f0^2+.03);
mh=real(ifftn(ifftshift(mf.*H)));
ShowSections(mh);
WriteMRC(mh,s.pixA,'emd_flat.mrc');
