% ShiftRotAlignMaps
% Rotate about the Z-axis and shift along the Z-axis to align maps of
% Cn-symmetric molecules
% m1c, m2c, msk1 (prelim. mask for m1c)

m1cr=rsRotateImage(msk1.*m1c,34);
n=size(m2c);
figure;
mysubplot(121);
imags(squeeze(sum(m1cr,3)));
msk=fuzzymask(n,3,n*.45,.1);
m2cm=msk.*m2c;

mysubplot(122);
P=Simplex('init',[34 1]);
%%
for i=1:100
    m2cr=rsRotateImage(m2cm,P(1)); %
    fsh=FourierShift(n,[0 0 P(2)]); % shift only in z
    m2crs=real(ifftn(fsh.*fftn(m2cr)));
    imags(squeeze(sum(m2crs,3)));
    drawnow;
    diff=m2crs-m1cr;
    err=diff(:)'*diff(:);
    P=Simplex(err);
    disp(P);
end;

