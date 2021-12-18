% P=ShiftRotAlignMaps(m1Ref,m2,P0);
% Rotate about the Z-axis and shift along the Z-axis to align maps of
% Cn-symmetric molecules
% m1Ref, m2, and P0 = [rot zshift]
% rot is in degrees, shift in pixels.

% Here's how I have used it to align I, apo and ReO maps
NISMapLoader5; % load the I map into m5, define msk1
m0=m5; % this is the NISMapLoader5 map, iodide.
m05=rsRotateImage(Crop(m5,144),35); % put it into standard position
m1Ref=Crop(m05,96).*msk1;
% NISMapLoader6;
% P0=[33.96 2.73] % 
% m2=m5;
% % or
% NISMapLoader7;
% P0=[67 2.3];
% m2=m5;
% % I've hardwired these values into NISMapLoader 6 and NISMapLoader7.
%
% exec
n=size(m1Ref);
figure;
mysubplot(121);
imags(squeeze(sum(m1Ref,3)));
msk=fuzzymask(n,3,n*.45,.1);
m2cm=msk.*Crop(m2,n);

%
mysubplot(122);
P=Simplex('init',P0);
%
for i=1:100
    m2cr=rsRotateImage(m2cm,P(1)); %
    fsh=FourierShift(n,[0 0 P(2)]); % shift only in z
    m2crs=real(ifftn(fsh.*fftn(m2cr)));
    imags(squeeze(sum(m2crs,3)));
    drawnow;
    diff=m2crs-m1Ref;
    err=diff(:)'*diff(:);
    P=Simplex(err);
    disp(P);
end;
P=Simplex('centroid')


