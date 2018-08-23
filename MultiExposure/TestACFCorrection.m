% TestACFCorrection.m
% uses m, acf


nd=size(ccf,1);
f=zeros(nd*nd,2);


% set up the linear least-squares
f(:,1)=ccf(:);
dcf=200*fuzzymask(nd,2,7.1,.01);
f(:,2)=dcf(:);
a=zeros(2,2);
for i=1:2
    for j=1:2
        a(i,j)=f(:,i)'*f(:,j);
    end;
end;



n=size(m,1);

cc=fftshift(real(ifftn(fftn(m(:,:,1)).*conj(fftn(m(:,:,2))))));

% cc=Crop(ccf+dcf*.1,n);

ccs=Crop(cc,ndis);




z=ccs(:);
y=f'*z;

coeffs=a\y;


cc2s=ccs-coeffs(1)*ccf;

% imacs(cc2s);
% 
% return
% 
% 
% 
% innprod=ccs(:)'*ccf(:);
% scale=innprod/(ccf(:)'*ccf(:));
% cc2s=ccs-scale*ccf;

%%
cc2=cc-Crop(ccf,n)*coeffs(1);
cc2s=Crop(cc2,ndis);
subplot(2,2,1);
imacs(ccs);
subplot(2,2,3);
plot(ccs(:,ndis/2+1));
subplot(2,2,2);
imacs(cc2s);
subplot(2,2,4);
plot(cc2s(:,ndis/2+1));