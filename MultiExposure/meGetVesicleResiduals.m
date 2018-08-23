function vesCorr=meGetVesicleResiduals(mi,mv)
% function vesCorr=meGetVesicleResiduals(mi,mv)
% Given the mi structure and the merged image mv with vesicles subtracted,
% compute a correction that can be subtracted from mv.  It's best if mv has
% already been highpass filtered, e.g. at .003 A^-1

% % test code
% fHP=.002; % Highpass, in A^-1
% cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1')
% % cd('/Users/fred/EMWork/Hideki/140620')
% d=dir('Info');  % entry 13 is sq10_001
% id=12+9;
% d(id).name
% mi=ReadMiFile(['Info/' d(id).name]);
% mi.basePath=AddSlash(pwd);
% m=meReadMergedImage(mi);
% n=size(m);
% disp('Make model vesicles');
% v=meMakeModelVesicles(mi,n);
% ds=mi.imageSize(1)/n(1);
% pixA=mi.pixA*ds;
% mv=GaussHP(m-v,fHP*pixA);
% displayOn=1;
% % end of test code

tRise=15; % Sharpness of radial mask, in A
rRise=.2; % sharpness of Fourier filter in radial direction
rRes=4;   % resolution of radial function

% display parameters
displayOn=0;
fLP=.2;
fLP2=.2;
dispExp=.5;

if displayOn
    figure(1);
    SetGrayscale;
    subplot(2,4,1);
    mvf=GaussFilt(mv,.05);
    disp('Start');
    tic
end;

n=size(mv);
ds=mi.imageSize(1)/n(1);
pixA=mi.pixA*ds;

vesCorr=zeros(n,'single');
goodVesInds=find(all(mi.vesicle.ok(:,1:3),2));

% nves=numel(mi.vesicle.x);

for ind=goodVesInds'
    %%
    %     Pick the working (extracted) image size
    r=mi.vesicle.r(ind,1)/ds;
    nw=128*ceil((1.2*pixA*r+80)/(pixA*64));  % size of extracted box
    ct=nw/2+1;
    pt=[mi.vesicle.x(ind) mi.vesicle.y(ind)]/ds+1;
    ipt=round(pt);
    %
    %     Extract the vesicle image
    mx=ExtractImage(mv,ipt,nw);
    mx=mx-median(mx(:));
    
    
    %     Convert to polar coordinates
    mp=ToPolar(mx);  % mp is nw/2 x 2*nw in size (radial x angular)
%     Take FT along angular dimension
    fmp=circshift(fft(mp'),[nw 0]);  % fftshift of angular fft
    
    %     Define the limits in the polar FT to be considered
    fw=max(2,round((r*pixA/80)));  % more cycles when r>200
    mw=(50+max(r*pixA-170,0)*.3)/pixA; % Radius half-width, in pixels
    rmsk=fuzzymask(nw/2,1,mw,tRise/pixA,r);  % 1D radius mask
    fR=pixA/rRes;
    if fR<1  % small pixels, we smooth the radial function
        rfilt=circshift(fuzzymask(nw/2,1,fR*nw/4,fR*rRise*nw/4),[nw/4 0]);
        fmp=ifft(fft(fmp').*repmat(rfilt,1,2*nw))';
    end;
    %      compute the correction
    fcorr=zeros(2*nw,nw/2,'single');  % angular x radial
%     Copy the masked radial function
    fcorr(nw-fw+1:nw+fw+1,:)=fmp(nw-fw+1:nw+fw+1,:)...
        .*repmat(rmsk',2*fw+1,1);
%     ifft of angular function
    cp=real(ifft(circshift(fcorr,[-nw 0])))';
    cr=ToRect(cp);
%     Add our small panel into the big vesCorr image.
    vesCorr=vesCorr+ExtractImage(cr,ipt,n,1);
    
    if displayOn
%           Show the overall micrograph with the position marked
        subplot(2,4,1);
        imacs(mvf);
        hold on;
        plot(pt(1),pt(2),'y+');
        hold off;
        title(mi.baseFilename,'interpreter','none');
        
%         Show the correction
        subplot(2,4,2);
        imacs(cr);
        %     plot(sect(cr),'b.-');

        %     Draw the extracted image with cross-hairs at the center
        mxf=GaussFilt(mx,fLP2);
        subplot(2,4,3);
        imacs(mxf);
        hold on;
        plot([ct ct nan 1 nw],[1 nw nan ct ct],'y-');
        hold off;
        title(['Vesicle ' num2str(ind) ';  r = ' num2str(round(r))]);
        
        %     Show it in polar coords
        subplot(2,4,7);
        imacs(mp');
        fmp(nw+1,1:2)=0;
        title('Polar coords');        
        
        %     Show the corrected image
        subplot(2,4,4);
        %     imacs(mx-cr);
        imacs(GaussFilt(mx-cr,fLP2));
        hold on
        plot([ct ct nan 1 nw],[1 nw nan ct ct],'y-');
        hold off;
        title('Corrected');

        % Show the corrected image in polar coords
        subplot(2,4,8);
        imacs(mp'-cp');
        title('Residual');

        %     Show the polar power spectrum
        smp=abs(fmp).^2;
        subplot(2,4,5);
        imacs(-10:10,0:nw/2-1,smp(nw-9:nw+11,:).^dispExp);
        ri=round(r);
        rmn=ri-mw+1;  % minimum radius to blank
        rmx=min(nw/2,ri+mw+1);  % max radius to blank
        hold on;
        plot(0,r,'y+');
        plot([-fw -fw fw fw -fw],[rmn rmx rmx rmn rmn],'r-');
        hold off;
        title(['Polar spectrum, dr = ' num2str(mw) ';  df = ' num2str(fw)]);
        
        subplot(2,4,6);
        spc=abs(fcorr).^2;
        imacs(-10:10,0:nw/2-1,spc(nw-9:nw+11,:).^dispExp);
        
        drawnow;
    end;
    
end;
%
%%
if displayOn
    toc;
%     subplot(231);
%     imacs(GaussFilt(m-v,.1));
%     subplot(233);
%     imacs(vesCorr);
%     subplot(232);    
    figure(2);
    clf;
    SetGrayscale;
    imacs(GaussFilt(mv,.3));
    pause(1);
    imacs(GaussFilt(mv-vesCorr,.3));

end;