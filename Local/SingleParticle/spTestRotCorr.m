% spTestRotCorr


%%  Test a single image pair.
img=GaussFilt(randn(64,64),.2)+0*randn(64,64);
ref=img+randn(64,64);
img(:,33)=img(:,33)+1;  % horizontal line
img=grotate(img,pi/180*10);
pimg=gridToPolar(img);
pref=gridToPolar(ref);
cc0=spRotCorr(pimg,pref);

nt=size(pimgs,2);
ccx=zeros(nt,1);
nt
tic
for t=1:nt
    q=grotate(ref,(t-1)*2*pi/nt);
    ccx(t)=q(:)'*img(:);
end;
toc
plot([cc0 ccx]);

%%

% Test lots of images
nim=100;
img=GaussFilt(randn(64,64),.2);
imgs=single(zeros(64,64,nim));
refs=single(zeros(64,64,nim));
for i=1:nim
    imgs(:,:,i)=grotate(img,(i-1)*0.2)+2*randn(64,64);
    refs(:,:,i)=grotate(img,(i-1)*.01)+randn(64,64);
end;
disp('gridToPolar');
tic
pimgs=gridToPolar(imgs);
prefs=gridToPolar(refs);
toc
disp('spRotCorr');
tic
ccs=spRotCorr(pimgs,prefs);
toc
imovie(ccs);

