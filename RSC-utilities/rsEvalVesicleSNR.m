% rsCountVesicles2
% Determine vesicle sizes and amplitudes
% from the Info directory of an experiment.
% Compare vesicle amplitudes^2 with noise spectra.

%workingDir='/Users/fred/EMWork/Hideki/160909p/KvLipo121_2w11v3m1/';
% workingDir='/Users/fred/EMWork/Hideki/160909p/KvLipo121_2w10v3';
% workingDir='/Users/fred/EMWork/Hideki/170808p/SimpleVes_raFit/';
% workingDir='/Users/fred/EMWork/Hideki/170808p/SimpleVes_w10/';
%workingDir='/Users/fred/EMWork/Hideki/180226/Kv_1sel_VesPW/';
% workingDir='/Users/fred/EMWork/Hideki/SNR/180226/Kv_1sel/';

%cd(workingDir);

startEntry=1
endEntry=5000

countEff=0.8;
showHist=0;
spectrumCorrection=1;
% spectrumCorrection=1/(1.53^2);
stride=1;
txtSize=12;
requireParticles=0;
d=dir('Info/');
imageFileSuffix='ms.mrc';
vesicleFileSuffix='mvs.mrc';
nPicks=0;
names=f2FindInfoFiles;
endEntry=min(numel(names),endEntry);
nEntries=endEntry-startEntry+1;

ds=4;  % assume we're working with downsampled images
nmi=0;
miPicks=zeros(nEntries,1);
miDef=zeros(nEntries,1);
vesR=cell(nEntries,1);
vesS=cell(nEntries,1);
vesOk=cell(nEntries,1);
miAmps=zeros(nEntries,1);
fs=(.005:.005:.1)';
nfs=numel(fs);
miSpecs=zeros(nEntries,nfs);
miShots=zeros(nEntries,1);
miDoses=zeros(nEntries,1);
miDefoci=zeros(nEntries,1);

nMod=100;
for nmi=1:nEntries
    name=names{startEntry+nmi-1};
    
    mi=ReadMiFile(name);
    if isfield(mi,'vesicle') && isfield(mi.vesicle,'x') && numel(mi.vesicle.x)>0 ...
            && (~requireParticles || (isfield(mi.particle,'picks') ...
            && numel(mi.particle.picks)>0))
        % we have vesicles with particles
        
        miDef(nmi)=mi.ctf(1).defocus;
        if numel(mi.vesicle.x)>1
            rawAmps=real(mi.vesicle.s(:,1));
            okVes=all(mi.vesicle.ok(:,1:3),2) & ~isnan(rawAmps) & rawAmps>0;
            if sum(okVes)>0
                miAmps(nmi)=median(rawAmps(okVes));
            end;
        end;
        %     fs=(.001:.001:.1)';
        [spec, shot, noiseOk]=meEvalNoiseModel(fs,mi);
        if noiseOk
            miSpecs(nmi,:)=spec';
            miShots(nmi)=shot(end);
        else
            disp([num2str(nmi) ' ' name ' -- No noise model.']);
        end;
        miDoses(nmi)=mi.doses(1);
        miDefoci(nmi)=mi.ctf(1).defocus;
        %     semilogy([spec shot]);
        if mod(nmi,nMod)==0
            disp([num2str(nmi) '  ' name]);
            nMod=max(nMod,10^(floor(log10(nmi))));
        end;
    else
        disp([num2str(nmi+startEntry-1) ' ' name ' -- No vesicles.']);
    end;
end;
%%
miShots=spectrumCorrection*miShots;
miSpecs=spectrumCorrection*miSpecs;
% return
%%
% function MakeDisplay
finds=[4 10 20];
nBins=60;
    legends=cell(numel(finds),1);
    for i=1:numel(finds)
        legends{i}=num2str(fs(finds(i)));
    end;


medDose=median(miDoses);
medShot=median(miShots(1:nmi));
miShots=min(max(miShots,medShot/10),medShot*10);
miMax=3*max(median(miSpecs(1:nmi,finds(1))),medShot);
miMin=.3*min(median(miSpecs(1:nmi,finds(end))),medShot);
oks=miSpecs(1:nmi,finds(end))>miMin & miSpecs(1:nmi,finds(1))<miMax;
[h,bins]=hist(log10(miSpecs(oks,finds)),nBins);
h1=0*h;
lmShot=log10(median(miShots(oks)));
ind=max(1,min(nBins,interp1(bins,1:nBins,lmShot,'nearest')));
h1(ind)=max(h(:));
if showHist
    subplot(211);
    % h1=hist(log10(miShots(oks)),bins);
    bar(bins,h);
    hold on;
    bar(bins,h1,'k');
    hold off;
    textLength=numel(workingDir);
    titleChars=40;
    xlabel('log LF spectra at given freqs');
    titleText=[workingDir(max(1,textLength-titleChars):end) '  (' num2str(medDose,2) ' e/�2)'];
    title(titleText,'interpreter','none');
    legend([legends;{'shot'}]);
    subplot(223);
else
    titleText=[(pwd) '  (' num2str(medDose,2) ' e/A^2)'];
%    mysubplot(311)
mysubplot(3,1,1,0,.05,0,0);
end;
miOkDoses=miDoses(oks);
nzDoses=miOkDoses>0;
miEstShots=0*miOkDoses;
miEstShots(nzDoses)=1./(countEff*miOkDoses(nzDoses)');
miEstShots(~nzDoses)=0.1;  % default value
inds=startEntry:endEntry;
semilogy(inds(oks),[miSpecs(oks,finds) miShots(oks) miEstShots]);
ylabel('shot, LF spectra');
legend([legends;{'shot'; 'estShot'}],'location','southeast');
title(num2str(median([miSpecs(oks,finds) miShots(oks)]),3))
if showHist
    subplot(224);
else
%    mysubplot(312);
mysubplot(3,1,2,0,.05,0,0);

end;
ampScl=1e3;
med=median(miAmps(oks));
plot(inds(oks),(miAmps(oks)*ampScl))
ymax=2*(med*ampScl);
axis([inds(1) inf 0 ymax]);
ylabel(['Vesicle amp x ' sprintf('%1.0e',ampScl)]);
text(endEntry*.02,ymax*.02,['Median vesicle amplitude ' num2str((med),3)],'verticalalignment','bottom');
title(titleText,'interpreter','none');

% mysubplot(nr,nc,ind,xs,ys,xo,yo)
mysubplot(3,1,3,0,.05,0,0);
plot(inds(oks),miDefoci(oks));
ylabel('Defocus, \mum');

drawnow
% end

%
% for i=1:1
%     msName=[mi.procPath mi.baseFilename imageFileSuffix];
%     [msName,ok]=CheckForImageOrZTiff(msName);
%     mvsName=[mi.procPath mi.baseFilename vesicleFileSuffix];
%     [mvsName,ok2]=CheckForImageOrZTiff(mvsName);
%     if ~(exist(msName,'file') && exist(mvsName,'file'     ))
%         disp('no file');
%         continue;
%     end;
%     m=ReadEMFile(msName);
%     n=mi.imageSize(1)/ds;
%     m=Downsample(m,n);
%     mv=ReadEMFile(mvsName);
%     mv=Downsample(mv,n);
%     v=m-mv;
%     mysubplot(221);
%     imags(mv);
%     mysubplot(222);
%     imags(v);
%     marksPresent=false;
%     nv=0;
%     if isfield(mi.vesicle,'x') && numel(mi.vesicle.x)>0 ...
%             && (~requireParticles || (isfield(mi.particle,'picks') ...
%             && numel(mi.particle.picks)>0))
%         % we have vesicles with particles
%         nv=numel(mi.vesicle.x);
%         disp([name '  ' num2str(nv) ' vesicles.']);
%         xs=double(mi.vesicle.x/ds+1);
%         ys=mi.vesicle.y/ds+1;
%         yls=double(ys+mi.vesicle.r(:,1)/ds);
%         effAmps=zeros(nv,1)+NaN;
%         rs=zeros(nv,1);
%         as=zeros(nv,1);
%         rs=zeros(nv,1);
%         oks=false(nv,1);
%         for j=1:nv
%             if any(mi.vesicle.ok(j,:)==0)
%                 continue;
%             end;
%             oks(j)=true;
%             effAmps(j)=1e4*EstimateImageAmplitude(mi,j,ds);
%             as(j)=mi.vesicle.s(j,1);
%             rs(j)=mi.vesicle.r(j,1)*mi.pixA;
%             disp(num2str([j effAmps(j) 1e4*as(j,1) rs(j)],3));
%         end;
%     end;
%     if nv>0
%         as1=as;
%         as1(as==0)=NaN;
%         subplot(2,2,3);
%         plot(rs,as,'bo');
%         subplot(2,2,4);
%         hist([as effAmps],50);
%         pause;
%     end;
% end;
%
%
% %%
%
% %         disp(' ');
%
%
% %         if size(mi.particle.picks,1)>0
% %             flags=mi.particle.picks(:,3);
% %             num=sum(flags>=16 & flags<48);
% %             nPicks=nPicks+num;
% %             miPicks(nmi)=num;
% %         end;
% %         vesR{nmi}=mi.vesicle.r;
% %         vesS{nmi}=mi.vesicle.s;
% %         vesOk{nmi}=mi.vesicle.ok & ~(requireParticles & size(mi.particle.picks,1)<1);
% %         if mod(nmi,stride)==0 || i>nEntries-3  % print out every 'stride' entries
% %             disp([num2str(nmi) ' ' d(i).name '  ' num2str([num nPicks])]);
% %             disp(abs(vesR{nmi}(1,:)));
% %         end;
% %
% % end; % if numel(mi.vesicle.x)
% % disp(' ');
% % if b=='Q'
% %     break;
% % end;
% % end; % for i
% % % miPicks=miPicks(1:nmi);
% % % miDef=miDef(1:nmi);
% % % subplot(313)
% % % hist(miPicks,100);
% % % xlabel('Particles per micrograph');
% % % ylabel('Frequency');
% % % subplot(312)
% % % plot(miPicks);
% % % ylabel('Particles per micrograph');
% % % xlabel('Micrograph');
% % % subplot(311);
% % % plot(miDef);
% % % ylabel('Defocus, \mum');
% % % xlabel('Micrograph');
% % disp('Done.');
% %
% %
% % return
%
% % %%
% % % Estimate vesicle amplitude from image integral
% % k=47;
% % ds=4;
% % ves=meMakeModelVesicles(mi,mi.imageSize(1)/ds,k,0,0);
% % imags(ves);
% % sv=sum(-ves(:))
% % a0=sum(mi.vesicleModel);
% % a1=sqrt(2*pi)*mi.vesicle.extraSD*(mi.vesicle.extraPeaks*0+1);
% % as=[a0 a1];
% % s=squeeze(mi.vesicle.s(k,1,:));
% % estA=4*pi*(mi.vesicle.r(k,1))^2*(as*s)/ds^2
% %
% % effVesAmp=ds^2*estA/(4*pi*(mi.vesicle.r(k,1)^2)*a0)
% % effVesAmp0=ds^2*sv/(4*pi*(mi.vesicle.r(k,1)^2)*a0)
%
%
% function effAmp=EstimateImageAmplitude(mi,k,ds)
% % Estimate image amplitude from normalized image
% ves=meMakeModelVesicles(mi,mi.imageSize(1)/ds,k,0,0);
% sv=sum(ves(:));
% mi1=mi;
% mi1.vesicle.s(k,:,:)=0;
% mi1.vesicle.s(k,1,1)=1;
% ves1=meMakeModelVesicles(mi1,mi.imageSize(1)/ds,k,0,0);
% sv1=sum(ves1(:));
% effAmp=sv/sv1;
% end
%
% % %% Compute coherence function
% % vpc0=meMakeModelVesicles(mi,960,0,1,1);
% % v1=mv+vpc0;
% % mi1=mi;
% % mi1.vesicle.s(:,:,2:3)=0; % delete the extra peaks
% % vp0=meMakeModelVesicles(mi1,960,28,0,0);
% % fv=fftshift(fftn(v));
% % fv0=fftshift(fftn(vp0));
% % %%
% % k=1e6;
% % H=(fv.*fv0)./(abs(fv0).^2+k);
% % figure(2);
% % SetComplex;
% % imacx(GaussFilt(H,.05),.3);
