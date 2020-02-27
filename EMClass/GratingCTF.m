% GratingCTF
fs=18;
sMax=.5;
s=(0:.0002:sMax)';
lambda=EWavelength(300);
Cs=0;
% Cs=2.7
B=0;
alpha=0;
% alpha=.05;
defocus=.25;

clf;
figure(5);
h=gcf;
h.Color='k';

[ctf,chi]=ContrastTransfer(s, lambda, defocus, Cs, B, alpha);
ys=ctf;
ys(:,2)=chi;

txt={'Contrast transfer' ; '\chi, radians'};
colorString={'r-' ; 'b-'};
for i=1:2
    subplot(2,1,i);
ha=gca;
ha.XColor='w';
ha.YColor='w';
ha.FontSize=fs;
ha.Color='k';
hold on;
plot(s,0*s,colorString{i},'linewidth',1.5);
plot(s,ys(:,i),'w-','linewidth',1.5);
xlabel('Spatial frequency s');
ylabel(txt{i});
text(0.9,sMax/10,'CTF = sin(\chi)','color',[1 1 .5],'fontsize',fs);
switch i
    case 1
        title(['Defocus ' num2str(defocus) ' \mum'],'color',[1 1 .8]);
    otherwise
        title('CTF = sin(\chi)','color',[1 1 .8]);
end;
hold off;
end;

