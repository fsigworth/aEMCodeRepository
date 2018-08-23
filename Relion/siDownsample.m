function sid=siDownsample(si,ds)
%function sid=siDownsample(si,ds)
% downsample an si structure by ds, changing pixA and all dimensions based
% on it; and cropping the ctfs and pwfs.
n=size(si.ctfs,1);
n1=n/ds;
sid=si;
sid.ctfs=Crop(si.ctfs,n1,1);
sid.ctf1s=Crop(si.ctf1s,n1,1);
sid.pwfs=Crop(si.pwfs,n1,1);
sid.yClick=si.yClick/ds;
sid.rVesicle=si.rVesicle/ds;
sid.pixA=si.pixA*ds;
sid.mbnOffset=si.mbnOffset/ds;

