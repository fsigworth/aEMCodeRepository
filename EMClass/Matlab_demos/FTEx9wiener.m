% FTEx9-wiener
% Compare phase-flipping with a Wiener filter for CTF-correction.

load smiley  % loads sm, fsm, c
snrfactor=100;
defocus=2;
B=40;
pixA=5;

noiseamp=20;

[np np1]=size(sm);
c=CTF(np,pixA,.025,defocus,2,B,0.1);  % 2 um defocus
c1=c(np/2+1:np,np/2+1);
 
fn=sft(noiseamp*randn(np));

fsmn=sft(sm);
smn=sm;

figure(1);
set(gcf,'color',[0.9 0.9 0.9]);
SetComplex;
nr=2; nc=4;

% subplot(nr,nc,1);
% imacs(smn);
% title([num2str(defocus) ' \mum defocus']);

subplot(nr,nc,1+nc);
% imacx(fsmn,0.4);
imacx(c);
title('CTF');

fsmc=fsmn.*c+fn;
smc=real(isft(fsmc));
subplot(nr,nc,1);
imacs(-smc);
title([num2str(defocus) ' \mum defocus']);


fpsmc=fsmc.*sign(c);
psmc=real(isft(fpsmc));
subplot(nr,nc,2);
imacs(psmc);
title('Phase flipped');

subplot(nr,nc,nc+2);
imacx(fpsmc,0.5);

% Get the average noise power
nvar=sum(sum(abs(fn).^2))/np^2;
% ivar=ToPolar(abs(fpsmc).^2,np/2,np*2,1,33,33);
ivar=ToPolar(abs(sft(sm)).^2,np/2,np*2,1,33,33);
tvar=sum(ivar,2)/(2*np);
svar=max(0,tvar)/nvar;
svar=svar*snrfactor;

w1=abs(c1).*svar./((c1.^2).*svar+1);  % Wiener weights as a function of radial frequency
w=ToRect(w1,np,1,[np/2+1 np/2+1]);
subplot(nr,nc,nc+3);
imacs(w);
title('Filter coefficients');

wfsmc=real(isft(fpsmc.*w));
subplot(nr,nc,3);
imacs(wfsmc);
title('Wiener filtered');

subplot(nr,nc,4);
plot(abs(c1));
title('|CTF|');

subplot(nr,nc,nc+4);
plot(w1);
title('Filter coefficients');

