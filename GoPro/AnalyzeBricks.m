cd('/Users/fred/Movies/GoPro/2018-08-10/HERO6 Black 1')
img=imread('GOPR0023.JPG');
m=single(img(:,:,1));
md=Downsample(m,[384 512]);
me=mean(md(:));
% hypot(384,512) = 640
mc=Crop(md,720,0,me);
imags(mc);

mr=grotate(mc,.975); mr(:,363)=0; imags(mr);
drawnow;
mr=grotate(mc,.975); mr1=mr(:,363); plot(mr1);
x=1:720;
thr=90+20/720*x+20/(720^6)*x.^6; plot([mr1 thr']);
mr2=thr'-mr1; plot(mr2);
%%
mrt=mr2>0;
plot(mrt)

mrt(1:75)=0;
plot(mrt)
ncrp=640;
mrtc=Crop(mrt,ncrp);
plot(mrtc)

%% Find center points
m1=mrtc;
p1=0;
nc=0;
centers=[];
while any(m1)
  p1=find(m1>.5,1);
  m1(1:p1)=1;
  p2=find(m1<.5,1);
  m1(1:p2)=0;
  nc=nc+1;
  centers(nc)=(p1+p2)/2;
end;

plot(centers);


%% angles at center points
xs0=.26;
nLeft=find(xs0+centers<ncrp/2,1,'last');
xs=xs0+(-nLeft:nc-nLeft-1)*20.06/69;  % tangents
angs=atan(xs);
subplot(221);
plot(angs);
grid on;

subplot(222);
plot(centers);

subplot(223);
plot(angs,[centers', 230*angs'+ncrp/2]);

subplot(224);
% resid=centers-240*angs+180*angs.^3-200*angs.^5+50*angs.^7);
% resid=centers-234.5*angs+11*sin(3.1*angs)-3*sin(5*angs)+0*sin(9.3*angs);
resid=[centers' (320+212.3*angs+18.5*angs.^3)']; % error is +/- .5 over most range.
% resid=centers-235*angs;  % error is +/- 10
plot(angs,resid);
grid on


subplot(222); % compute the inverse, angles from pixels
resangs=angs-centers/216+0.05*((centers-320)/212).^3;
resangs=[angs' (-1.485+(centers/216-0.05*((centers-320)/212).^3)')];
plot(centers,resangs);
%%
% What is angular max?  

% mrtr=mrtc+flipud(circshift(mrtc,[12,0]));
% plot(mrtr)
% 
% mrts=mrtr>1.5;
% 
% plot(mrts);
% 
% diagonal pixel
mxPix=hypot(512,384)/2;
mxAng=mxPix/220*180/pi;  % about 83 degrees.








da=atan(20.03/69);  % base angle, in radians
