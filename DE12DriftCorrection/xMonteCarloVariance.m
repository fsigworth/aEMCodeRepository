% xMonteCarloVariance
% Conclusion is that precision of location goes as SNR for SNR (amp/SD)
% > 0.25, but as sqrt(SNR) below that corner.  precision (val/SD) is down
% to 1/2 at SNR=.07 (val/SD).

n0=1e6;
n=1e6;
nas=200;
a0s=zeros(nas,1);
vars=zeros(nas,1);
for i=1:nas
    a0s(i)=.001*10^((i-1)*5/nas);
    %     n=round(n0/sqrt(a0s(i)));
    a=a0s(i)+randn(n,1);
    b=randn(n,1);
    angs=atan2(b,a)/(2*pi);
    s1=mean(angs);
    s2=mean(angs.^2);
    vars(i)=s2;
    if mod(i,10)==0
        hist(angs,1000);
        title(i);
        drawnow;
    end;
end;
%%
effst=1./(1./vars-12);
% loglog(a0s,vars);
fit=2*pi*a0s.*(1+(1.2./a0s).^2.5)./(1+(1.3./a0s).^2).*exp(-.3./a0s);
% the above is a good fit down to about SNR=.3, and then drops off rapidly.
% fit=2*pi*a0s.*exp(-.4./a0s);
% the above is a conservative fit that starts to drop at SNR=1.
loglog(a0s,sqrt(1./effst),a0s,a0s*(2*pi),a0s,sqrt(a0s)*(pi),a0s, fit);
axis([-inf inf 1e-4 inf]);
ylabel('effective sd of position');
xlabel('snr (amp / noise sd)');
