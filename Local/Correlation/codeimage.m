% codeimage.m
% Try cross-correlation of two noisy images

load modelVesicle.mat;
numShifts=100;
n=size(modelVesicle);
sigmaN1=1;  % noise standard deviation
sigmaN2=.1;
noisyVesicle1=modelVesicle+sigmaN1*randn(n); %adding noise for this trial
noisyVesicle2=modelVesicle+sigmaN2*randn(n);
sumtrial=zeros(numShifts,1);
dxs=zeros(numShifts,1);

for i=1:numShifts
    dx=i-ceil((numShifts+1)/2);
    dxs(i)=dx;
    dy=0;
    shiftedVesicle=circshift(noisyVesicle2,[dx dy]); %shifting the noisy model vesicle by dx= 50 circularly
    %     a=noisyVesicle.*shiftedVesicle;  % multiplying each element of the noisy vesicle array and the shifted noisy vesicle array
    %     b=cumsum(a(:)); %finding and making an array the cumulative sum of each product of the two arrays
    %     sum=b(end); %taking the last element of the cumulative sum.
    sum=noisyVesicle1(:)'*shiftedVesicle(:);  % simple inner product
    sumtrial(i)=(1/(n(1)*n(2)))*sum;  %multiplying this sum by 1/N*M
    %     fprintf('value with shift %d : %d\n', dx, sumtrial(i)) %printing out the sum for trial i.
end
[pkval,pki]=max(sumtrial);  % Get the peak value and the index
pkx=dxs(pki);  % peak location, in pixels, away from the center.

% sumtrial=transpose(sumtrial);
figure(1);
colormap gray(256)
subplot(2,2,1);
imagesc(noisyVesicle1);

subplot(2,2,2);
imagesc(noisyVesicle2);

subplot(2,2,3);
plot(dxs,sumtrial,'r',pkx,pkval,'ko') %plots the 100 trials and the peak.
xlabel('X-displacement');
ylabel('correlation');
title(['Peak location : ' num2str(pkx)]);

