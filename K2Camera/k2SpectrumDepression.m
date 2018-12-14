% function spd=k2SpectrumDepression(fs)
% d=2.2;  % modified from Li et al, JSB '13, where it is 2.58
% spectrum should be 1-b*spd, where b=-.0194+.0295*dr (Fitted from Table 1
% values) and dr is dose rate in e/pixel-s
% 
% d=2.2;
%     spd=sinc(fs*d).^2;
% 
% a=125.5;
% b=.17;
% d=2.2;
% 
% fs=(1:n/2)/n;
% 
% S=a*(1-b*(sinc(fs*d)).^2)';


rates=[1.02 2.02 4.05 8.16 9.73 12.87 14.9 20.16];
bs=[.0048 .043 .0991 .222 .271 .363 .423 .570];

plot(rates,bs);
p=polyfit(rates,bs,1);
p
bc=polyval(p,rates);
plot(rates, [bs' bc' rates'*.03-.02]);
