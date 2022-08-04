% FTEx8spectrum
% Show the power spectrum of a random object operated on by the ctf.

% load smiley  % loads sm, fsm, c
% [np np1]=size(sm);
np=128;

defocus=2;  % microns
c=CTF(np,5,.025,defocus,2,100,0.1);  % 2 um defocus

smn=32*randn(np);
fsmn=sft(smn);

figure(1);
SetComplex;
set(gcf,'color',[0.5 0.5 0.5]);
set(gcf,'menubar','none');

nr=2; nc=4;

subplot(nr,nc,1);
imags(smn);
title('Random object');
axis off


subplot(nr,nc,1+nc);
imacsx(fsmn,0.4);
title('FT of object');
axis off


subplot(nr,nc,2);
imags(isft(c));
title('Point-spread');
axis off


subplot(nr,nc,2+nc);
imacsx(c);
title('CTF');
axis off


fsmc=fsmn.*c;
smc=real(isft(fsmc));

psmc=real(isft(fsmc));
subplot(nr,nc,3);
imags(psmc);
title('Image');
axis off


subplot(nr,nc,nc+3);
imacsx(fsmc,0.5);
title('FT of image');
axis off


spsmc=abs(fsmc).^2;
ispsmc=real(isft(spsmc));
ispsmc(np/2+1,np/2+1)=0;

subplot(nr,nc,4);
imags(sign(ispsmc).*(abs(ispsmc).^0.5));
title('ACF');
axis off


subplot(nr,nc,nc+4);
imacs(sqrt(spsmc));
title('Power spectrum');
axis off


