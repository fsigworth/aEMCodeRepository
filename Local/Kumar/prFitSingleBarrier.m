% prFitSingleBarrier
% Fit single barrier to prestin data
% mV, icl, isc
 cd('/Users/fred/Documents/papers/Kumar prestin')
%  load Fig3Data.mat
 load Fig3DataSD2.mat
 
for k=1:2
scnMode=k-1;
if scnMode
    isc=prISCN-ctISCN;
    sd=sqrt(sdCtISCN.^2+sdPrISCN.^2);
    P=[ .5  1   0   0   .75  4  1];
    msk=[1  1   0    0   0  1  0]; % scn
    y=isc*1000;
    ey=sd*1000;
else
    icl=prICl-ctICl;
    sd=sqrt(sdCtICl.^2+sdPrICl.^2);
    P=[ .5 .5   0  0  .75  1  1];
    msk=[1  1   0   0  0  0  0]; % cl
    %   d1  g1  d2  g2 q pr1 pr2
    y=icl*1000;
    ey=sd*1000;
end;

subplot(2,1,scnMode+1);

for j=1:2
P=Simplex('init',P,.1,msk);
for i=1:1500
    delta1=P(1);
    delta2=P(3);
    g1=P(2);
    g2=P(4);
    q=P(5)/25;
    pr1=P(6);
    pr2=P(7);
    I=g1*(pr1*exp(mV*q*delta1)-exp(-mV*q*(1-delta1)))+g2*(pr2*exp(mV*q*(1-delta2))-exp(-mV*q*delta2));
    diff=(y-I);
    diff(1:2)=0;
    diff(end-1:end)=0;
    err=diff'*diff;
    err=err+1e6*sum(P<0)+1e6*(P(1)>1)+1e6*(P(3)>1);
    if mod(i,10)==0
    plot(mV,y,'.',mV,I,'-','markersize',20);
    drawnow;
    end;
    P=Simplex(err);
end;
end;
%%
%     plot(mV,y,'.',mV,I,'-','markersize',20);
    plot(mV,I,'k-');
hold on;
errorbar(mV,y,ey,'.','markersize',20);
hold off;
grid on
hold off;
% title(num2str(P));
xlabel('Membrane potential, mV');
ylabel('Current density pA/pF');
SubplotLabel(0.07,0.85,char(64+k));

P

G1=g1*25/P(5)*1e-3
G2=g2*25/P(5)*1e-3

end;  % for k