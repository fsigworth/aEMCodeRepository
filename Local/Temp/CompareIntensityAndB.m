% CompareIntensityAndB.m
% See how fitted B values depend on various parameters

% save Ps1Quadrants Ps Intens 
load Ps1Quadrants % gets Ps and Intens
load Mis
load IimgsF12.mat % gets Iimgs

nim=numel(Mis);
nq=size(Ps,2);
npi=size(Iimgs,1);
nxq=floor(npi/nq);

quadMeans=zeros(nim,nq^2);
quadSDs=zeros(nim,nq^2);

% analyze Iimgs
for i=1:nim
    for qi=1:nq^2 % quadrant index
        qix=mod(qi-1,nq)+1;
        qiy=floor((qi-1)/nq)+1;
        im=Iimgs((qix-1)*nxq+1:qix*nxq,(qiy-1)*nxq+1:qiy*nxq,i);
        quadMeans(i,qi)=mean(im(:));
        quadSDs(i,qi)=std(im(:));
    end;
end;



%%

sdThresh=.025;  % Maximum std for acceptable quadrants
bThresh=110;
iThresh=.7;
%

nFs=zeros(nim,1);
for i=1:nim
    nFs(i)=diff(Mis{i}.frameSets(1,:))+1;
end;

% Pick up defocus, B and intensity for each quadrant.
Ds=zeros(nim,nq,nq);
Bs=zeros(nim,nq,nq);
Is=zeros(nim,nq,nq);
for i=1:nim
    for j=1:nq
        for k=1:nq
            Ds(i,j,k)=Ps(i,j,k).defocus;
            Bs(i,j,k)=Ps(i,j,k).B;
            Is(i,j,k)=Intens(i,j,k)/nFs(i);
        end;
    end;
end;
I2=reshape(Is,nim,nq^2);
% I2=reshape(Intens/16.5,nim,nq^2);
B2=reshape(Bs,nim,nq^2);
D2=reshape(Ds,nim,nq^2);

clf;

Imax=Percentile(I2,.97);

ok=B2<bThresh&quadSDs<sdThresh&I2>iThresh*Imax;
B20=B2(ok);
I20=I2(ok);
D20=D2(ok);
nB=sum(ok(:));

B20n=B2;
B20n(~ok)=NaN;
I20n=I2;
I20n(~ok)=NaN;
D20n=D2;
D20n(~ok)=NaN;

subplot(221);
plot(I20n/Imax,B20n,'.');
axis([0.7 1.1 0 bThresh]);
xlabel('Intensity');
ylabel('B');
title('200 images, 3x3 quadrants');

% Do a lin fit
F=zeros(nB,2);
F(:,1)=ones(nB,1);
F(:,2)=I20(:);
a=LinLeastSquares(F,B20);
hold on;
plot(I20/Imax,F*a,'k-');
legend(num2str(a(2)*Imax),'-0.1 intens -> B+12');
hold off;


subplot(223);
plot(I2/Imax,D2,'.');
xlabel('Intensity');
ylabel('Defocus');
axis([0.7 1.1 -inf inf]);

% Do a lin fit
% a=LinLeastSquares(F,D20);
slope=.29;
a=[-slope*Imax slope]';
b=1.4:.2:3;
y1=repmat(F*a,1,9);
b1=repmat(b,nB,1);
y=kron(y1,b);
hold on;
% plot(I20,(I20-Imax)*.28+1.95,'k-');
plot(I20/Imax,y1+b1,'k-');
legend(num2str(a(2)*Imax),'-0.1 intens -> d-75nm','location','northwest');
hold off;
% decrease in intensity to .75 <-> .2 um defocus change.
% .075um change in defocus per 10% intensity change.

% There's no correlation between defocus and B
subplot(222);
plot(quadSDs,B2,'.',quadSDs*0+sdThresh,B2,'k-');
axis([-inf inf 0 bThresh]);
xlabel('Quadrant intensity SD');
ylabel('B');

% plot(D20n,B20n,'.');
% xlabel('Defocus');
% ylabel('B');
% axis([1.2 3 0 100]);
% 
subplot(2,2,4);
plot(quadMeans,quadMeans*0+sdThresh,'k-',...
     quadMeans,quadSDs,'.');
xlabel('Intensity');
ylabel('Quadrant intensity SD');

