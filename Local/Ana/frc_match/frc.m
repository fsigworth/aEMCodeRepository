function coeff=frc(I1,I2,rings)
nRings=size(rings,3);
coeff=zeros(1,nRings);
I1=fftshift(fft2(I1));
I2=fftshift(fft2(I2));
for i=1:nRings,
    ring=squeeze(rings(:,:,i));
    tmp1=I1.*ring;
    tmp1=tmp1(:);
    tmp2=I2.*ring;
    tmp2=tmp2(:);
    coeff(i)=(real(tmp1'*tmp2))/(norm(tmp1,2)*norm(tmp2,2));
end