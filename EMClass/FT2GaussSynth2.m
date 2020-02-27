% FT2GaussSynth.m
% Example of building up a 2d image from Fourier components.
%
vlSetDarkGraphics(20);
clf;
vlSet1080Figure(1,.75,1120);
mkSize=10;
nx=256;
ndis=5;
ndis2=16;
ctr=nx/2+1;
dctr=ndis2/2+1;
xLim=2;
useQuadrant=1;
makeVideo=1;

cd('/Users/fred/Documents/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.4FT2D/Figs/')
vName='FT2GaussSynth';
if makeVideo
    v=VideoWriter(vName);
    v.FrameRate=2;
    open(v);
    disp(['Making movie file: ' vName '.avi']);
end;
x=(-xLim:2*xLim/nx:xLim-2*xLim/nx)';
disX=-ndis2/2:ndis2/2-1;
numel(x)
R=RadiusNorm(nx)*xLim*2;
m=exp(-pi*R.^2);
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
% disSquare=Crop(fm,ndis2);
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
        fVal=disQuadrant(j,i);
        if abs(fVal)<1e-2
            continue;
        end;
        
        hold on;
        i1=i-1;
        j1=j-1;
        if useQuadrant
            for ip=-1:2:1
                for jp=-1:2:1
                    fTerm1(ctr+i1*ip,ctr+j1*jp)=fVal;
                end;
            end;
        else
            fTerm1(ctr+i1,ctr+j1)=fVal;
            fTerm1(ctr-i1,ctr-j1)=conj(fVal);
        end;
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
        term1=real(fftshift(ifftn(ifftshift(fTerm1))))*(nx/2/xLim)^2;
        subplot(236);
        imaga(x,x,term1*512+128);
        
        str=sprintf('(%d,%d)  %2.4f',i-1,j-1,fVal);
        title(str,'color','w');
        
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
        
        if makeVideo
            f=getframe(gcf);
            writeVideo(v,f);
        end;
    end;
end;
if makeVideo
    close(v);
end;



