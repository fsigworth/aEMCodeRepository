% BFactorPlot.m
% estimate how many more images we need for higher resolution
% k*exp(-B.^2/kDiv)

kDiv=2; % 2 for power, 4 for amplitude
res=(4.5:-.5:2.5)';
f=1./res;
Bs=50:50:100;
nb=numel(Bs);
a=repmat(f,1,nb);

for i=1:nb
a(:,i)=exp(-Bs(i)*f(:).^2/2).*res(:)./res(1);
end;

plot(res,a);
disp([res a]);