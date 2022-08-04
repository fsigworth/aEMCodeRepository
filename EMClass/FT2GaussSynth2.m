% FT2GaussSynth.m
% Example of building up a 2d image from Fourier components.
%
vlSetDarkGraphics(20);
clf;
vlSet1080Figure(1,.75,1120);
mkSize=10;
nx=256;
ndis=5;
ndis2=11;
ctr=nx/2+1;
% dctr=ndis2/2+1;
xLim=2;
useQuadrant=0 ;
disScale=1000;
makeVideo=1;


% cd('/Users/fred/Documents/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.4FT2D/Figs/')
cd('/Users/fred/Documents/Documents - Katz/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.4FT2D/Figs/')
vName=['FT2GaussSynth2' char(98-useQuadrant)]; % 2a means quadrant, 2b means whole disc
%  of Fourier space.
if makeVideo
    v=VideoWriter(vName,'MPEG-4');
    v.FrameRate=2;
    open(v);
    disp(['Making movie file: ' vName '.mpr']);
end;
% x=(-xLim:2*xLim/nx:xLim-2*xLim/nx)';
x=(-xLim:2*xLim/(nx-1):xLim)';  % include endpoints
disX=-(ndis2-1)/2:(ndis2-1)/2;

numel(x)
R=RadiusNorm(nx)*xLim*2;
m=exp(-pi*R.^2);
% Display the original Gaussian
subplot(232);
imags(x,x,m);

% Make half-sized plots on the left
subplot(231);
ha=gca;
pos1=ha.Position;
pos2=pos1;
pos2(1)=pos1(1)+pos1(3)*.8;
pos2(3)=pos1(3)*.4;
ha.Position=pos2;
plot(sect(m,2),x);
ha.YTickLabel={};
ha.XTickLabel={};


msh=ifftshift(m);
fm=fftshift(real(fftn(msh)))/(nx/2/xLim)^2;

disQuadrant=fm(ctr:ctr+ndis,ctr:ctr+ndis);
disSquare=Crop(fm,ndis2);
xdis=(0:1+1/ndis:ndis)';
subplot(233);
if useQuadrant

    imags(xdis,xdis,disQuadrant);
else
    imags(disX,disX,disSquare);
    % colormap jet
    % imagesc(xdis,xdis,(disQuadrant));
    % axis xy
    % colorbar;
end;

% subplot(223);
% imags(disQuadrant>1e-6)

termSum=zeros(nx,nx);
iOld=NaN;
jOld=NaN;
dctr=ndis2/2+1;

ndis1=6;
for j=1:ndis1
    for i=1:ndis1
        subplot(233);
        hold on;
        plot(iOld,jOld,'bo','markersize',mkSize);
        hold off;
        
        fTerm1=zeros(nx,nx);
        fTerm1d=zeros(nx,nx); % displayed ft
        fVal=disQuadrant(j,i);
        if abs(fVal)<1e-2
            continue;
        end;
        
        hold on;
        i1=i-1;
        j1=j-1;
%         put in points and Friedels
            fTerm1d(ctr+i1,ctr+j1)=fVal;
            fTerm1d(ctr-i1,ctr-j1)=conj(fVal);
            
            
            for ip=-1:2:1
                for jp=-1:2:1
%                     sel=2*ip+jp+5; % selection index runs 1..4
                    fTerm1(ctr+i1*ip,ctr+j1*jp)=fVal; % actual points we fill in.
                    plot(i1*ip,j1*jp,'go','markersize',mkSize);
                end;
            end;

%             fTerm1=fTerm1d;
%         end;
        plot(i1,j1,'ro','markersize',mkSize);
        if ~useQuadrant
            plot(-i1,-j1,'ro','markersize',mkSize);
        end;
        hold off;
        
        iOld=i1;
        jOld=j1;
        
        %         if all([i j]==1)
        %             fVal=fVal/4;
        %         elseif any([i j]==1)
        %             fVal=fVal/2;
        %         end;
        % Compute the IFT of the 
        term1d=real(fftshift(ifftn(ifftshift(fTerm1d))))*(nx/2/xLim)^2;
        subplot(236);
        imaga(x,x,term1d*disScale+128);
        
        str=sprintf('(%d,%d)  %2.4f',i-1,j-1,fVal);
        title(str,'color','w');
        
%         Compute the FT of the transposed points too, and use that for
%         reconstruction
        term1=real(fftshift(ifftn(ifftshift(fTerm1))))*(nx/2/xLim)^2;
        termSum=termSum+term1;
        subplot(235);
        %         plot(sum(termSum,2)/nx);
        imaga(x,x,termSum*260);
        
        % Make half-sized plots on the left
        subplot(234);
        ha=gca;
        pos1=ha.Position;
        pos2=pos1;
        pos2(1)=pos1(1)+pos1(3)*.8;
        pos2(3)=pos1(3)*.4;
        ha.Position=pos2;
        plot(sect(termSum,2),x);
        ha.YTickLabel={};
        ha.XTickLabel={};
        drawnow;
        
        if makeVideo
            f=getframe(gcf);
            writeVideo(v,f);
        else
            pause;
        end;
    end;
end;
        hold on;
        plot(iOld,jOld,'bo','markersize',mkSize);
        hold off;
if makeVideo
    close(v);
    disp(['Movie ' vName ' written.'])
    disp(['In: ' pwd])
end;



