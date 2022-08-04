% FTEx7-CTF
%
% Demonstrate CTF effect on an image, and demonstrate phase-flipping

load smiley  % loads sm, fsm, c

defocus=3;

c=CTF(64,5,.025,defocus,2,100,0.1);  % 2 um defocus, 5 Å pixels

figure(1);
SetComplex;
nr=2; nc=4;

subplot(nr,nc,1);
imags(sm);

subplot(nr,nc,1+nc);
imacsx(fsm,0.4);

subplot(nr,nc,2);
imags(isft(c));

subplot(nr,nc,2+nc);
imacsx(c);

fsmc=fsm.*c;
smc=real(isft(fsmc));
subplot(nr,nc,3);
imacs(-smc);
title('CTF, inverted');

subplot(nr,nc,3+nc);
imacsx(fsmc,0.5);

fpsmc=fsmc.*sign(c);
psmc=real(isft(fpsmc));
subplot(nr,nc,4);
imacs(psmc);
title('Phase flipped');

subplot(nr,nc,nc+4);
imacsx(fpsmc,0.5);
