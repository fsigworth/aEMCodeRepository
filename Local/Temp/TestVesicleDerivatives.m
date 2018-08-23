% TestVesicleDerivatives.m

dr=1;


for i=1:3;
    r=100+(i-2)*dr;
    v(:,:,i)=VesicleFromModel(n,r,[1 1]);
end;
dv=v(:,:,3)-v(:,:,1);
gv=v(:,:,1)+v(:,:,3)-2*v(:,:,2);
imacs(GaussFilt(gv,.05));
% imacs(v(:,:,2))
