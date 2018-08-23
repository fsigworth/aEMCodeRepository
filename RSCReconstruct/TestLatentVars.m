% TestLatentVars.m
% Look at latent variables when stopping reEMStep26.
% Assume nt=5, nAlpha=3

figure(4);
q=reshape(diffIR,5,5,6,338);
q1=q(3,3,[2 5],:);
q1=squeeze(q1)';
plot(q1)
ylabel('diffIR');
subplot(222);
imags(q(:,:,2,158))
plot(q1)
r=rawLogProbs.priors(:,:,:,20);

subplot(223);
r=reshape(r,5,5,6,338);
r1=r(3,3,[2 5],:);
r1=squeeze(r1);
plot(r1')
ylabel('prior');
subplot(224);
imags(r(:,:,2,158));
xlabel('prior at rso, 158');

figure(5);
subplot(221);
p=reshape(ccIR,5,5,6,338);
p1=squeeze(p(3,3,[2 5],:));
plot(p1')
ylabel('ccIR');

subplot(222);
imags(p(:,:,2,171));
xlabel('ccIR at rso,171');
imags(p(:,:,2,158));

subplot(223);
ctRefImgs=zeros(1024,nRV);
for i=1:nRV, ctRefImgs(pix(:),i)=ctRefs(:,i);end; ctRefImgs=reshape(ctRefImgs,32,32,nRV);
imags(ctRefImgs(:,:,158));

rImgs=zeros(1024,150);
for i=1:150, rImgs(pix(:),i)=transImgs(:,i);end; rImgs=reshape(rImgs,32,32,5,5,6);
subplot(224);
imags(rImgs(:,:,25+13))

figure(8)
imags(M)
colormap jet
imacs(ccIR')
[val,i,j]=max2d(ccIR')
val =
    4.8468
i =
   164
j =
    38
imacs(diffIR');
