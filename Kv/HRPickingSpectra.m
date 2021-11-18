% HRPickingSpectra
% Make figures of power spectra and SNR for Kv TM maps (Figs. 6,7)
% First we calculate with no CTF at various B values. (Fig. 8)
% then at B=80 for various CTFs.

load HRPicking/Figs/TM_Micelle_maps.mat % mMap0 tMapsh compMap1sh micDens1
n2=180;
tmOnly=1;
hfCutoff=1/(6*1.06); % effect of downsampling by 3 of original micrograph
nDis=n2;

%% old version: without and with micelles.
% if tmOnly
%     map=Crop(tMapsh,n2);
% else
%     map=Crop(compMap1shx,n2);
% end;
% micMap=Crop(mMap0,n2);
% 
% 
% %%
% % See the total power in acf as a function of resolution.
% disp('Computing radial power spectra...');
% sp2=RadialPowerSpectrum(map);
% sp2(:,2)=RadialPowerSpectrum(map+micMap);
% disp('...upsampled spectra...')
% sp4=RadialPowerSpectrum(Crop(map,n4));
% sp4(:,2)=RadialPowerSpectrum(Crop(map+micMap,n4));
% disp('done.')

%% new version: without and with T1 domain (mapT, mapC)
    mapT=Crop(tMapsh,n2);
    mapC=Crop(compMap1shx,n2);
maps=mapT;
maps(:,:,:,2)=mapC;

% n2=size(tMapsh,1);
us=2;
hp1=.06;
n4=n2*us;
% ctr2=ceil((n2+1)/2);
pixA=1; % 1 angstrom voxels
freqs2=(0:n2/2-1)'/n2;
freqs4=(0:n4/2-1)'/n4;
n2t=find(freqs2>.35); % high frequency limit
n2z=find(freqs2>hp1,1); % find first frequency above .06
freqs2t=(0:n2t-1)'/n2;
n4z=find(freqs4>.06,1) % find first frequency above .06
n4t=find(freqs4>.35,1)
freqs4t=(0:n4t-1)'/n4;
cShifts=[0 0 0; 0 0 -21];

%%
% See the total power in acf as a function of resolution.
disp('Computing radial power spectra...');
sp2=RadialPowerSpectrum(mapT);
sp2(:,2)=RadialPowerSpectrum(mapC);
disp('...upsampled spectra...')
sp4=RadialPowerSpectrum(Crop(mapT,n4));
sp4(:,2)=RadialPowerSpectrum(Crop(mapC,n4));
disp('done.')


sp2=sp2(1:n2t,:); % truncate at n2t
sp20=sp2;
sp2(1:n2z,:)=0; % Force lf points to zero.
sp4=sp4(1:n4t,:); % truncate at n4t
sp40=sp4;
sp4(1:n4z,:)=0; % Force lf points to zero.




%% images plus spectra
figure(5);
labels={'TM' 'TM+T1'};
for i=1:2
    mysubplot(2,4,i*4-3);
    proj=squeeze(sum(circshift(maps(:,:,:,i),cShifts(i,:)),2));
    imags(rot90(Crop(proj,nDis),2));
    hold on;
    plot([20 120],[5 5],'w-','linewidth',2);
    text(70,7,'100 �','color','w','horizontalalignment','center','verticalalignment','bottom');
    text(5,n2-5,labels{i},'color','w','fontsize',16,'horizontalalignment','left','verticalalignment','top');
    axis off 
end;
subplot(1,2,2);
semilogy(freqs4t,sp40,'linewidth',1);
ylabel('Spectral density S_S');
xlabel('Spatial frequency, Å^{-1}');
grid on;
limits=[0 .3 1e-6 1e-1];
ys=limits(3:4)';
xs=[.06 1/6];
hold on;
plot([xs;xs],[ys ys],'-','color',[.5 .5 .3]);
hold off;
legend('TM','TM+T1','f_{low}','f_{high}');
axis(limits);


    
%%



%% Show the SNR in the absence of CTF

figs=[6 7];
scaleUp=1000;

yMax=[3.5 8]*scaleUp/100;
bVals=[0 50 100];
nBs=numel(bVals);
freqs4pin=repmat(freqs4t*2*pi,1,nBs);
x016=find(freqs4t>1/6,1);

labels={'B=0'};

    oldSNR=zeros(n4t,nBs);
oldSSN=zeros(n4t,nBs);


for i=1:size(sp4,2) % do for each spectrum type
    sps=sp40(:,i);
    for k=2:nBs % modify the spectra according to B
        B=bVals(k);
        labels{k,1}=['B=' num2str(B)];
        gaussB=exp(-freqs4t.^2*B/2);
        sps(:,k)=sps(:,1).*gaussB;
    end;
    % Compute the SNR as a function of frequency cutoff
%     cumPower=cumsum(sps.*freqs4pin,1).^2;
%     cumNoiseVar=cumsum(sps.*freqs4pin,1);
%     cumSNR=cumPower./cumNoiseVar;
    spsz=sps;
    spsz(1:n4z,:)=0;
    cumSNR=cumsum(spsz.*freqs4pin,1);
    cumSNR=cumsum(sps.*freqs4pin,1);
%     cumSpect=cumsum(sps,1);

    figure(figs(i));
    subplot(131);
    if i==1
        aLims1=[0 .4 0 .17];
        aLims2=[0 .4 0 6];
        
    else
        aLims1=[0 .4 0 .17];
        aLims2=[0 .4 35 55];
        
    end;
    
        proj=squeeze(sum(circshift(maps(:,:,:,i),cShifts(i,:)),2));

    imags(rot90(Crop(proj,nDis),2));
    axis off 

%     subplot(222);
%     semilogy(freqs4t,sps,'linewidth',1.5);
%     legend(labels,'Location','northeast');
%     ylabel('Mean spectral density');
%     grid on
% %     axis([0 .35 1e-8 1.2]);
%     xlabel('Spatial frequency, Å^{-1}')

    rings=2*pi*repmat(freqs4t,1,nBs);
    subplot(132);
    plot(freqs4t,sps.*rings*scaleUp,'linewidth',1.5);
    legend(labels,'Location','northeast');
    ylabel('S \times 2\pi f   (arbitrary units)');
    xlabel('Spatial frequency, Å^{-1}')
    grid on
     axis(aLims1);

    subplot(133);
    plot(freqs4t,cumSNR*scaleUp,'linewidth',1.5);
%     hold on;
%     set(gca,'colororderindex',1);
%     plot(freqs4t,oldSNR);
%     hold off;
    legend(labels,'location','northwest');
    ylabel('Relative SNR');
    grid on
%     axis([0 .35 0 yMax(i)]);
    xlabel('Spatial frequency, Å^{-1}')
%     axis(aLims2);

end;



%% --------------

%% Find signal as a function of defocus
% SNR at B=0 is given by sp1 x 2\pi r
% freqs4=(0:n4/2-1)'/n4;
% n4t=find(freqs4>.35,1);
% freqs4t=freqs4(1:n4t);
% find first frequency above .06
% n4t=52*us; % truncated frequencies
% sp40t=sp4(1:n4t,:);
% sp40t(1:8*us,:)=0; % zero before f=.05;
bVals=80;
% freqs4pi=freqs4t*2*pi;
figs=[8 9];
% k=1; % TM only
% defs=[.19 .31];
defs=[.31 1 2];
nDefs=numel(defs);
signals=zeros(n4t,nDefs);

ctPars.lambda=.0197;
ctPars.alpha=.07;
ctPars.Cs=2.7;
ctPars.B=bVals(1);
ctAlt=ctPars;
ctAlt.alpha=pi/2;
ctAlt.Cs=0;
ctAlt.defocus=0;

freqs4piz=freqs4pin(:,1);
freqs4piz(1:n4z)=0;

cumSNR=zeros(n4t,nDefs);
spectScls=[7000 6000];

for k=1:2 % TM, Comp maps
    figure(figs(k));
    subplot(141);
    
    proj=squeeze(sum(circshift(maps(:,:,:,k),cShifts(k,:)),2));
    imags(rot90(Crop(proj,nDis),2));
    axis off
    
    
    
    subplot(142);
    cla;
    hold on;
    % set(gca,'colororderindex',4);
    for i=1:nDefs
        %     if defs(i)~=0
        ctPars.defocus=defs(i);
        ct2=ContrastTransfer(freqs4t,ctPars).^2;
        %     else
        %         ct2=ContrastTransfer(freqs4t,ctAlt).^2;
        %     end;
        if i==1
            lw=1.5;
        else
            lw=1;
        end;
        if i<nDefs
            plot(freqs4t,ct2,'linewidth',lw);
        end;
            cumSNR(:,i)=cumsum(ct2.*sp40(:,k).*freqs4pin(:,1)*scaleUp);
        cumSNRz(:,i)=cumsum(ct2.*sp40(:,k).*freqs4piz*scaleUp);
        
    end;
    
    plot(freqs4t,sp40(:,k)*spectScls(k),'k--','linewidth',1);
    hold off;
    grid on;
    axis([0 inf 0 1])
    legTexts=cell(nDefs,1);
    for i=1:nDefs
        legTexts{i}=['d =' num2str(defs(i)) '\mum'];% ' S=' num2str(cumSNR(end,i),2)];
    end;
    legTexts2=legTexts;
    legTexts2{end}='Spectrum';
    legend(legTexts2);
    xlabel('Spatial frequency, Å^{-1}')
    ylabel('CTF^2 and Power Spectrum');
    
    subplot(143);
    % plot(freqs4t,signals(:,1));
    %  set(gca,'colororderindex',4);
    % hold on;
    plot(freqs4t,cumSNR,'linewidth',1);
    legend(legTexts,'location','southeast');
    % hold off;
    xlabel('Spatial frequency, Å^{-1}')
    ylabel('Relative SNR');
    grid on;
    
    subplot(144);
    % plot(freqs4t,signals(:,1));
    %  set(gca,'colororderindex',4);
    % hold on;
    plot(freqs4t,cumSNRz,'linewidth',1);
    if k==2
        hold on;
        plot(.16,cumSNRz(x016,2),'k+','markersize',10,'linewidth',1.5);
        hold off;
        legTexts{end+1}='test condition ';
    end;
    legend(legTexts,'location','southeast');
    
    % hold off;
    xlabel('Spatial frequency, Å^{-1}')
    ylabel('Relative SNR;  f > .06 only');
    grid on;
end;

%     legend(legTexts(1:end),'location','best');

% subplot(222);
% plot(defs,signals);

% %% Look at the individual radial spectra of the projections
% figure(5); clf;
% sp1s=RadialPowerSpectrum(Crop(projs,n0,1),1);
% %
% sp1sd=sp1s.*repmat(freqs2',1,npP);
% plot(freqs2,sp1sd);
% % sort angles by tilt
% [vals,inds]=sort(abs(angs(:,2)-90));
% for i=1:4
%     lower=floor((i-1)*npP/4+1);
%     upper=floor(i*npP/4)
%     subplot(2,2,i);
%     plot(freqs2,sp1sd(:,lower));
% end;
% %
% for i=1:npP
%     plot(freqs2,sp1sd(:,i));
%     title(vals(i));
%     drawnow;
% end;
% 
