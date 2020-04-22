% GratingCTFPlots: plot CTF and \chi functions
vlSetDarkGraphics;
bkgColor=[1 1 1]*44/256;
fs=18;
sMax=.5;
s=(0:.0002:sMax)';
lambda=EWavelength(300);
Cs=0;
Cs=2.7;
B=0;
alpha=0.5;
alpha=.05;
defocus=.25;

clf;
figure(5);
h=gcf;
h.Color=bkgColor;

[ctf,chi]=ContrastTransfer(s, lambda, defocus, Cs, B, alpha);
ys=ctf;
ys(:,2)=chi;

txt={'Contrast transfer' ; '\chi, radians'};
colorString={'w--' ; 'w--'};
for i=1:2
    subplot(2,1,i);
ha=gca;
ha.XColor='w';
ha.YColor='w';
ha.FontSize=fs;
ha.Color=bkgColor;
hold on;
plot(s,0*s,colorString{i},'linewidth',1);
plot(s,ys(:,i),'w-','linewidth',1.5);
xlabel('Spatial frequency s');
ylabel(txt{i});
text(0.9,sMax/10,'CTF = sin(\chi)','color',[1 1 .5],'fontsize',fs);
switch i
    case 1
    if Cs==0
        title(['Defocus ' num2str(defocus) ' \mum;  Amplitude contrast ' num2str(alpha)],'color',[1 1 .8]);
    else
        title(['Defocus ' num2str(defocus) ' \mum;  Amplitude contrast ' num2str(alpha) ...
         '; C_s=' num2str(Cs) 'mm' ],'color',[1 1 .8]);
    end;
    otherwise
        title('CTF = sin(\chi)','color',[1 1 .8]);
end;
hold off;
end;

