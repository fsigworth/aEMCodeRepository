function vesCorr=meGetVesicleResidual(mi,mv)
% Given the mi structure and the merged image mv with vesicles subtracted,
% compute a correction that can be subtracted from mv.  It's best if mv has
% already been highpass filtered, e.g. at .003 A^-1

% test code
fHP=.003; % Highpass, in A^-1
cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1')
% cd('/Users/fred/EMWork/Hideki/140620')
d=dir('Info');  % entry 13 is sq10_001
id=12+9;
d(id).name
mi=ReadMiFile(['Info/' d(id).name]);
mi.basePath=AddSlash(pwd);
m=meReadMergedImage(mi);
n=size(m);
v=meMakeModelVesicles(mi,n);
mv=GaussHP(m-v,fHP*pixA);




tRise=15; % Sharpness of radial mask, in A

fLP=.2;
fLP2=.4;
dispExp=.5;

ds=mi.imageSize(1)/n(1);
pixA=mi.pixA*ds;

if displayOn
    figure(1);
    SetGrayscale;
    subplot(2,3,1);
    mvf=GaussFilt(mv,.05);
end;

vesCorr=zeros(n,'single');
nves=numel(mi.vesicle.x);
for ind=1:nves
    %%
    %     Pick the working (extracted) image size
    r=mi.vesicle.r(ind,1)/ds;
    nw=128*ceil(r/45);  % size of extracted box
    ct=nw/2+1;
    pt=[mi.vesicle.x(ind) mi.vesicle.y(ind)]/ds+1;
    ipt=round(pt);
    %
    %     Extract the vesicle image
    mx=ExtractImage(mv,ipt,nw);
    mx=mx-median(mx(:));
    
    
    %     Convert to polar coordinates
    mp=ToPolar(mx);  % mp is nw/2 x 2*nw in size (radial x angular)
    fmp=circshift(fft(mp'),[nw 0]);  % fftshift
    
    
    %     Determine the frequency range to be considered
    fw=max(2,round((r/25)));
    %     Make a mask in the radial variable
    mw=(50+max(r*pixA-170,0)*.3)/pixA;
    msk=fuzzymask(nw/2,1,mw,tRise/pixA,r);
    %      compute the correction
    fcorr=fmp.*0;
    fcorr(nw-fw+1:nw+fw+1,:)=fmp(nw-fw+1:nw+fw+1,:)...
        .*repmat(msk',2*fw+1,1);
    cp=real(ifft(circshift(fcorr,[-nw 0])))';
    cr=ToRect(cp);
    vesCorr=vesCorr+ExtractImage(cr,ipt,n,1);
    
    if displayOn
        subplot(2,3,1);
        imacs(mvf);
        hold on;
        plot(pt(1),pt(2),'y+');
        hold off;
        title(mi.baseFilename,'interpreter','none');
        
        %     Draw the extracted image with cross-hairs at the center
        mxf=GaussFilt(mx,fLP);
        subplot(2,3,2);
        imacs(mxf);
        hold on;
        plot([ct ct nan 1 nw],[1 nw nan ct ct],'y-');
        hold off;
        title(['Vesicle ' num2str(ind) ';  r = ' num2str(round(r))]);
        
        %     Show the power spectrum
        subplot(2,3,3);
        imacs(mp');
        fmp(nw+1,1:2)=0;
        title('Polar coords');
        
        %     Show the polar power spectrum
        smp=abs(fmp).^2;
        subplot(2,3,4);
        imacs(-10:10,0:nw/2-1,smp(nw-9:nw+11,:).^dispExp);
        ri=round(r);
        rmn=ri-mw+1;  % minimum radius to blank
        rmx=min(nw/2,ri+mw+1);  % max radius to blank
        hold on;
        plot(0,r,'y+');
        plot([-fw -fw fw fw -fw],[rmn rmx rmx rmn rmn],'r-');
        hold off;
        title(['Polar spectrum, dr = ' num2str(mw) ';  df = ' num2str(fw)]);
        %     Show the corrected image
        subplot(2,3,5);
        %     imacs(mx-cr);
        imacs(GaussFilt(mx-cr,fLP2));
        hold on
        plot([ct ct nan 1 nw],[1 nw nan ct ct],'y-');
        hold off;
        title('Corrected');
        subplot(2,3,6);
        imacs(cr);
        %     plot(sect(cr),'b.-');
        %     imacs(mp'-cp');
        title('Residual');
        drawnow;
    end;
    
end;
%

if displayOn
    subplot(231);
    imacs(GaussFilt(m-v,.1));
    subplot(232);
    imacs(GaussFilt(m-v-vesCorr,.1));
    subplot(233);
    imacs(vesCorr);
end;