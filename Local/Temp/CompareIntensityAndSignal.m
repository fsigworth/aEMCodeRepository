% CompareIntensityAndSignal
ds=32;      % Downsampling for intensity map
fc0=.0012;  % Filter for intensity map
iv0=2;
rmin=80;  % smallest radius we consider at all.
rmax=350;
rsel=[120 170];  % selected radius range
rBins=rmin:10:rmax;
onlyVesWithParticles=0;

imgOffset=0;

maxImg=10;

cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1');
miDir=dir('Info/');
ivals=[];
svals=[];
imeds=[];
smeds=[];
selRs=[];
allRs=[];
partAmps=[];
partAbsAmps=[];
medAbsAmps=[];
Iimgs=zeros(120,120,1);

figure(1);
SetGrayscale;

for fi=1:maxImg-imgOffset
    ind=fi+12+imgOffset;
    %     [fi ind]
    mi=ReadMiFile(['Info/' miDir(ind).name]);
    m0=ReadMRC([mi.imagePath mi.imageFilenames{1}]);
    md=BinImage(m0,ds);
    me=mean(md(:));
    msk=meGetMask(mi,mi.imageSize/ds);
    
    subplot(2,3,1);
    imacs(md);
    title(miDir(ind).name,'interpreter','none');
    
    intens=LaplacianFilter(md-me,msk,fc0*ds,100)+me;
    intensf=intens/size(mi.frameShifts{1},1);  % intensity per frame
    Iimgs(:,:,fi)=intensf;
%     subplot(2,3,2);
%     imacs(intens);
%     drawnow;
    % end;
    % save IimgsF12.mat Iimgs
    % return
    % for fi=1:0
    nv=numel(mi.vesicle.x);
    if nv>0
        sh=max(abs(mi.vesicle.s(:,2:end)),[],2);  % abs(higher order s terms)
        medSh=median(sh(sh>0));  % median of all nonzero values
        rh=max(abs(mi.vesicle.r(:,2:end)),[],2);  % abs(higher order radius terms)
        medRh=median(rh(rh>0));
        goodRh=rh<3*medRh & mi.vesicle.r(:,1)>rmin/mi.pixA; % small asymmetries, nonzero r
        goodRv=mi.vesicle.r(:,1)>rsel(1)/mi.pixA & mi.vesicle.r(:,1)<rsel(2)/mi.pixA;
        goodS=sh<3*medSh;
        goodOk=all(mi.vesicle.ok,2);
        good=   goodOk & goodS & goodRh & goodRv;
        goodSOk=goodOk & goodS & goodRh;  % not discriminated according to r value
        
        if isfield(mi,'particle') && isfield(mi.particle,'picks')...
                && numel(mi.particle.picks)>0
            goodParts=mi.particle.picks(:,3)>=16;
            amps=mi.particle.picks(goodParts,5);
            np=numel(amps);
            partVes=mi.particle.picks(goodParts,4);  % particles' vesicles
            vesicleBins=hist(partVes,1:nv)';
            hasParticle=vesicleBins>0;
            partAmps(np,fi)=0;
            partAmps(1:np,fi)=amps;
            absa=amps.*mi.vesicle.s(partVes,1);
            partAbsAmps(np,fi)=0;
            partAbsAmps(1:np,fi)=absa;
            partAbsAmps(np:end,fi)=NaN;
            subplot(2,3,2);
            hist(partAbsAmps(:),1e-4:1e-4:5e-3);
            xlabel('Particle amplitudes');
            medAbsAmps(fi,1)=median(absa);
        else
            hasParticle=false(nv,1);
            medAbsAmps(fi,1)=NaN;
            partAbsAmps(:,fi)=NaN;
        end;
        hasParticle=hasParticle | ~onlyVesWithParticles;
        good=good&hasParticle;
        goodSOk=goodSOk&hasParticle;
        if sum(good)>0
            rSelVals=mi.vesicle.r(good,1);
            rVals=mi.vesicle.r(goodSOk,1);
            xs=round((mi.vesicle.x(good)-1)/ds+1);
            ys=round((mi.vesicle.y(good)-1)/ds+1);
            xs=max(1,min(xs,size(intensf,1)));
            ys=max(1,min(ys,size(intensf,2)));
            ivals1=intensf(sub2ind(size(intensf),xs,ys));
            imed1=median(ivals1);
            svals1=mi.vesicle.s(good,1);
            smed1=median(svals1);
            n1=numel(ivals1);
            ivals(n1,fi)=0;
            ivals(1:n1,fi)=ivals1-iv0;
            ivals(n1+1:end,fi)=NaN;
            svals(n1,fi)=0;
            svals(1:n1,fi)=svals1;
            imeds(fi,1)=imed1;
            smeds(fi,1)=smed1;
            selRs(n1,fi)=0;
            selRs(1:n1,fi)=rSelVals;
            selRs(n1+1:end,fi)=NaN;
            n2=numel(rVals);
            allRs(n2,fi)=0;
            allRs(1:n2,fi)=rVals;
            allRs(n2:end,fi)=NaN;
            
            subplot(2,3,3);
            plot(ivals+iv0,svals,'.');
            
            subplot(2,3,4);
            plot(xs,ys,'k.');
            
            subplot(2,3,5);
            theRs=allRs(:);
            theRs(theRs==0)=NaN;
            theRs=theRs(~isnan(theRs));
            hist(theRs(:)*mi.pixA,rBins);
            axis([rmin rmax 0 inf]);
            xlabel('Vesicle radius, Å');
            
            subplot(2,3,6);
            sMin=.0025;
            hist(svals(~isnan(svals)&svals>sMin),sMin:.0002:.015);
            axis([sMin .015 0 inf]);
            xlabel('Vesicle amplitude');
            
            %         plot(selRs*mi.pixA,svals,'.');  % plot amplitude vs radius
            %         axis([rsel -inf inf])
            %         xlabel('Radius');
            %         ylabel('Amplitude');
            
        else
            disp(['no good vesicles, ind= ' num2str(ind)]);
        end;
    end;
    drawnow;
end;
imeds(imeds==0)=NaN;
ivals(ivals==0)=NaN;
ivals0=ivals+iv0;
%%
imeds(imeds==0)=NaN;
smeds(smeds==0)=NaN;

subplot(2,3,4);
plot(ivals0,svals,'.');
hold on
xm0=1.48;
xsp=9.5e-3;
ys=xsp*(imeds-xm0);
plot(imeds,ys,'k-');
hold off;
axis([min(ivals0(:)) max(ivals0(:)) .003 .013]);
xlabel('1st exposure dose');
ylabel('Vesicle amplitude');
ni=numel(imeds);
strs=cell(1,0);
if ni<10
    for i=1:ni
        strs{i}=num2str(i+imgOffset);
    end;
    legend(strs,'location','southeast');
end;

subplot(2,3,3);
q=[smeds ((imeds-xm0)*xsp) medAbsAmps*5.5];
% q=[smeds (imeds*2e-3)];
plot(imgOffset+1:imgOffset+size(q,1),q,'.-');
if onlyVesWithParticles
    title('Only vesicles with Particles');
end;
legend('signal','intensity fit','particle fit','location','southwest');
xlabel('Image number');
ylabel('Amplitude');
