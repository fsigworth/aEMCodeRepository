% FTEx0_1DFT.m

rectMode=1;  % pick rect -> sinc or gauss -> gauss
fSize=12;
afSize=14;
a=1;  % scaling factor

figure(1);
clf;
set(gcf,'toolbar','none');
set(gcf,'menubar','none');

maxX=2.5+2*rectMode;
xStep=.01;
xs=(-maxX:xStep:maxX)';
zs=0*xs;
subplot(121);
set(gca,'fontsize',14);
subplot(122);
plot(xs,zs,'k--');
axis([-maxX maxX -.5 1.5]);
xlabel('u');
ylabel('G(u)');
grid on;
set(gca,'fontsize',afSize);

Gs=[];
Xs=[];
nY=0;
g=zs;
if rectMode
    g(abs(a*xs)<.5)=a;
    G=sin(pi*xs)./(pi*xs);
else
    g=a*exp(-pi*(xs*a).^2);
    G=exp(-pi*(xs).^2);
end;
f=0;
b=0;
startup=1;
while b~='q'
    c=cos(2*pi*f*xs);
    subplot(121);
    plot(xs,[3.5+g 2.2+c 2.2+0*c c.*g]);
    hold on;
    plot(xs,[3.5+zs 2.2+zs zs],'k--');
    hold off;
    axis([-maxX maxX -1.2 6]);
    xlabel('x');
    ylabel('                    g(x)');
    set(gca,'fontsize',afSize);
    
    text(0.6*maxX,3.7,'g(x)','fontsize',fSize);
    text(0.6*maxX,2.4,'cos(2\piux)','fontsize',fSize);
    if b~=0 % startup
        nY=nY+1;
        Xs(nY)=f;
        Gs(nY)=xStep*(c'*g); % value of the integral
        text(0.6*maxX,0.2,num2str(Gs(nY),4),'fontsize',fSize);
        subplot(122);
        plot(Xs,Gs,'ko',xs,zs,'k--');
        axis([-maxX maxX -.5 1.5]);
        xlabel('u');
        ylabel('G(u)');
        grid on;
        set(gca,'fontsize',afSize);
    end;
    [f,a,b]=Myginput(1);
end;
hold on;
plot(xs,G);
% subplot(121);  % show we're done: clear the left panel.
% cla;

