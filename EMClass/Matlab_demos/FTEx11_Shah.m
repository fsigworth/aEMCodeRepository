% FTEx11  Shaw reconstruct.
disp('FTEx1: Show the FT of the Shah function');
% hit a key in Command Window to add each additional term to reconstruct

figure(1); set(gcf,'menubar','none','color',[0.8 0.8 0.8]);
n=256;
per=64;
% x=fftshift(0:n-1)';
x=(0:n-1)'/per;
y=ones(n,1);


for nt=1:per-1
    yadd=2*cos(2*pi*nt*x);
    y=y+yadd;
    plot(x,[y y*0]);
    axis([0 inf -per/3 2*per]);
    title(nt);
pause;
end;



% 
% % subplot(1,2,1); plot(x/10,y);
% % title('Function');
% % axis([-inf inf -.5 2]);
% Y=fft(fftshift(y));
% % subplot(1,2,2); plot(x/10,fftshift(real(Y))/10)
% % axis([-inf inf -.5 2]);
% % title('FT');
% 
% figure(3); set(gcf,'menubar','none');
% % show the reconstruction
% y1=zeros(1,n);
% k=2*pi/n;
% fy=real(Y*2/n);
% fy1=fy;
% fy1(1)=fy1(1)/2;
% for i=0:nterms
%     dy=fy1(i+1)*cos(i*k*x);
%     y2=y1+dy;
% %     subplot(1,2,1); plot([y1' y2' dy'-0.5]);
% %     legend(['sum ' num2str(max(i-1,0))],['sum ' num2str(i)],...
% %         ['term ' num2str(i)]);
%     subplot(1,2,1);
%     plot(x*dx,[y2' dy'-0.5]);
%     hold on;
%     plot(x*dx, [dy'*0 dy'*0-0.5],'k--');
%     hold off;
%     
% %     legend(['sum ' num2str(i)]);
%     y1=y2;
%     axis([-inf inf -1 1.5]);
% %     hold on;
% %     plot(dy-0.5);
% %     hold off;
%     title(['term ',num2str(i)]);
%     subplot(1,2,2);
%     plot(x2*dx,fy(1:n/2)*200/n);
%     hold on;
%     plot(x2(1:i+1)*dx,fy(1:i+1)*200/n,'ko');
%     hold off;
%     title('Fourier coefficients');
%     pause
% end;
% disp('Done.');
% 
