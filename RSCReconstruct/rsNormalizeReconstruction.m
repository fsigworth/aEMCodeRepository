function map=rsNormalizeReconstruction(rvol,nvol,k)
n=rvol.n;
ks=3;
comp=gridMakePreComp(n,ks);

vol=gridMakeNullFT(n,3);
vol.PadFT=rvol.PadFT./(k+nvol.PadFT);

% Get the reconstructed volume.
map=gridRecoverRealImage(vol,comp);

