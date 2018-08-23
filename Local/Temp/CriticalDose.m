% CriticalDose.m
% Critical dose curve at 300kV from Grant & Grigorieff, eLife '15

a=.245;
b=-1.665;
c=2.81;
kVFactor=.75;

x=0:.01:.25;
f=x;

Nc=.75*(.245*f.^-1.665+2.81);
% Nc=k*(a*x.^b+c);

plot(x,Nc);