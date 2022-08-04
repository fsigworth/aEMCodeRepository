% FTEx10.IFT
% Make a movie of Fourier synthesis of a Gaussian function.
% The movie is written out as a series of jpeg images.
% fs Nov 2018
writeImgs=0;
figure(1);
nt=6;
np=200;
w=np/4;
xs=(-2*w:2*w-1)'/w;
y=exp(-pi*xs.^2);
fy=real(fft(fftshift(y)))/(2*w);

subplot(122);
plot(xs,y);

sy=fy(1)/2+0*xs;
plot(xs,[y sy]);
terms=zeros(np,nt);
base=zeros(1,nt);
for i=1:nt
   subplot(121);
   terms(:,i)=cos(i/2*pi*xs)*fy(i+1);
   base(i)=nt/2-.35*i+2*abs(fy(i+1));
   bases=repmat(base(1,1:i),np,1);
   plot(xs,bases,'k--');
   hold on;
   plot(xs, terms(:,1:i)+bases,'b-');
   for j=1:i
       text(xs(1)+0.1,bases(1,j),num2str(fy(j+1),3),'verticalalignment','bottom');
   end;
   hold off;
   axis([-2 2 .8 nt-2]);
   set(gca,'ytick',[]);
   sy=sy+terms(:,i);
   subplot(122);
   plot(xs,[y sy]);
   axis([-2 2 -.2 1.1]);
   if writeImgs
   figName=[num2str(i) '.jpg'];
   disp(figName);
   print(figName,'-djpeg');
   end;
   end;