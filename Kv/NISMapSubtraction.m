% NISMapSubtraction.m


showRadialAverage=1; % show just the radial average.

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Silvia')

mapName1='NIS_I_Na_resampled.mrc';
mapName1='MapsForSubtraction/cryosparc_P1_J55_010_volume_map_sharp(6)_copy.mrc';
mapName2='MapsForSubtraction/NIS_Na.mrc';
pdbName1='MapsForSubtraction/SR_7_1_21_dimer_copy.pdb';

[m1,s]=ReadMRC(mapName1);
m2=ReadMRC(mapName2);
p1=pdb2mat(pdbName1);

% locate the ions
[u,a,b]=unique(p1.element);
ptrsNa=find(b==4);
ptrsI=find(b==2);
iPtrs=[ptrsNa(2:3); ptrsI(1)];
nIons=numel(iPtrs);

%%
n1=size(m1,1);
% n=144;
n=96;
us=2; % upsampling for rotate
vs=2; % further upsampling for plotting
nu=us*n;
nv=vs*nu;

msk1=fuzzymask(n,3,0.45*n,.1*n);

m1c=Crop(m1,n);
m1cr=rsRotateImage(msk1.*m1c,34);
m1uc=Downsample(m1c,n*us);

m2c=rsRotateImage(Crop(m2,n),-2);
m2cr=rsRotateImage(msk1.*Crop(m2,n),32);
m2uc=Downsample(m2c,n*us);

bigFigure=2;

% ?***
m1uc=m2uc;
bigFigure=3;
% ****

%%
% Make synthetic maps
muSy=zeros([1 1 1]*nu,'single'); % all atom positions marked
muSi=zeros([1 1 1]*nu,'single'); % ion positions marked
ctru=nu/2+1;
ctr0=n1/2;
na=numel(p1.X);

% Make coords = coordinates of atoms in the 2x upsampled map.
c0=[p1.X;p1.Y;p1.Z];
coords=round(us*(c0/s.pixA-ctr0))+ctru;

ok=all((coords>=1 & coords<=nu),1);
for i=1:na
    if ok
        muSy(coords(1,i),coords(2,i),coords(3,i))=1;
    end;
end;
for i=1:3
    j=iPtrs(i);
    muSi(coords(1,j),coords(2,j),coords(3,j))=1;
end;


%%
muSy=SharpFilt(muSy,.33/us);
m1Sy=Downsample(muSy,n);

muSi=SharpFilt(muSi,.33/us);
m1Si=Downsample(muSi,n);

%%
figure(1);
% imags(sum(m4c-600*m4Sy,3));
imags(sum(m1uc+1000*muSi,3));

%%
% m2c=Crop(m2,n);
% figure(2);
% ShowSections(m2c);
%%
uScale=500;
nnInterp=0; % 1: faster ERotate3 with nearest neighbor
showSi=1;

% rotations from standard position:
% none, 45deg about z, y, x

if showRadialAverage
    angs=[0 0 0];
    range=6*us*vs;
    
else
    
    angs=[0    0    0;
        pi/4 0    0;
        0   pi/4  0;
        pi/2 pi/4 -pi/2];
    range=5*us*vs;
end;
xVals=(-range:range)*s.pixA/(us*vs);
yLims=[-.51 1.51];
nAngs=size(angs,1);
lines=zeros(2*range+1,3,nAngs,nIons); % (dens, dirs, ion)

for i=1:nAngs
    
    phiAng=34;
    phi=phiAng*pi/180+angs(i,1);
    theta=angs(i,2);
    psi=angs(i,3);
    
    ctru=nu/2+1;
    m1ucr=ERotate3(m1uc,[phi theta psi],[1 1 1]*ctru,nnInterp);
    % Further upsapling
    m1vcr=Downsample(m1ucr,nv);
    
    if showSi % show synthetic ions
        musr=ERotate3(muSi,[phi theta psi],[1 1 1]*ctru,nnInterp);
        mvsr=Downsample(musr,nv);
    end;
    %     show the projection with ion positions marked
    figure(1);
    imags(sum(m1vcr,3));
    title(num2str(angs(i,:)*180/pi))
    
    % Crearte the rotated coordinates p1r
    % mrot=RotMatrix2(-psiAng);
    mrot=EulerMatrix(phi,theta,psi);
    % mrot(3,3)=1;
    ctr0=n1/2;
    p1Coords=[p1.X;p1.Y;p1.Z]/s.pixA-ctr0;
    p1rCoords=mrot*p1Coords;
    p1r.X=p1rCoords(1,:);
    p1r.Y=p1rCoords(2,:);
    p1r.Z=p1rCoords(3,:);
    
    ctrv=nv/2+1;
    hold on;
    plot(p1r.X(iPtrs)*us*vs+ctrv,p1r.Y(iPtrs)*us*vs+ctrv,'bo');
    hold off;
    
    %
    %
    % mysubplot(122);
    % plot(p1r.X,p1r.Y,'.');
    % axis([-n/2 n/2 -n/2 n/2]);
    drawnow;
    
    % Make plots of density along x,y,z for each ion.
    %%
    figure(bigFigure);
    legtxt=cell(3*nAngs,1);
    for ia=1:nAngs
        for ix=1:3
            legtxt{ix+(ia-1)*3}=[char(double('w')+ix) num2str(ia)];
        end;
    end;
    ltp=1:3*nAngs;
    for j=1:nIons
        ptr=iPtrs(j);
        ionLabel=[p1.element{ptr} ' chain ' p1.chainID{ptr}];
        %     Get coords into the upsampled map.
        cds=round([p1r.X(ptr) p1r.Y(ptr) p1r.Z(ptr)]*us*vs+ctrv);
        
        if showRadialAverage
            mysubplot(2,nIons,j);
            imags(m1vcr(:,:,cds(3))); % show the section
            hold on;
            plot(cds(1),cds(2),'yo','markersize',20);
            hold off;
            axis off;
            title(ionLabel);
            %             nx=NextNiceNumber(2*range+5)
            nx=2*range+3;
            mLoc=ExtractVolumeInterp(m1vcr,cds,nx);
            [rMean,rMedian,yVals,rVals]=Radial3(mLoc,[]);
%             % Correct the r=0 point
%             rMean(1)=0*rMean(1)+1*mean(yVals{2}(1:6)); % 1/2 times r=1 points
%             rMean(2)=4/3*rMean(2)-1/3*mean(yVals{2}(1:6)); % remove the r=1 points
            
            sMean=[rMean(range+1:-1:1); rMean(2:range+1)];
            sMedian=[rMedian(range+1:-1:1); rMedian(2:range+1)];
            lines(:,1,1,j)=sMean;
            lines(:,2,1,j)=sMedian;
            legTxt{j}=ionLabel;
            %             radLine=zeros(2*range+1);
            %             radLine(range+1:end)=rMed(1:range);
            %             radLine(range:-1:1)=rMed(2:range-1);
            %
            mysubplot(2,nIons,nIons+j)
            plot(xVals,[sMean sMedian]);
            plot(xVals,sMedian,'linewidth',1);
            hold on; 
            plot(xVals,0*xVals,'k-');
            hold off;
            grid on;
            axis([-inf inf yLims]);
            ylabel('Map density');
            xlabel('Radius, Å');
        else
            
            rngs0=[{0} {0} {0}];
            for k=1:3 % x,y,z trace
                rngsh=rngs0;
                rngsh(k)={-range:range};
                lines(:,k,i,j)=squeeze(m1vcr(cds(1)+rngsh{1},cds(2)+rngsh{2},cds(3)+rngsh{3}));
            end;
            %         mysubplot(2,nIons,j);
            %         linesToPlot=reshape(lines,2*range+1,3*nAngs,nIons);
            %         plot( (-range:range)*s.pixA/(us*vs) , linesToPlot(:,ltp,j) ,'linewidth',1);
            %         axis([-inf inf min(lines(:)) max(lines(:))]);
            %         legend(legtxt);
            % Show the vicinity of the ion
            iAng=mod(i-1,2);
            kAng=floor((i-1)/2);
            col=2*(j-1)+iAng+1;
            rowA=1+kAng;
            rowB=3+kAng;
            %         disp([i j rowB col]);
            nr=4;
            nc=nIons*2;
            mysubplot(nr,nc,nc*(rowA-1)+col);
            curves=1:3;
            if i>1
                lines(:,5-i,i,j)=0; % zero out the redundant trace
            end;
            plot( xVals, lines(:,:,i,j) ,'linewidth',1);
            hold on;
            plot( xVals , 0*xVals, 'k-','linewidth',.5);
            hold off;
            axis([-inf inf min(lines(:)) max(lines(:))]);
            if i==1 && j==1
                legend('x','y','z');
            end;
            grid on;
            if i==1 % first angle
                title(ionLabel);
            end;
            
            mysubplot(nr,nc,nc*(rowB-1)+col);
            imags(m1vcr(:,:,cds(3)));
            hold on;
            plot(cds(1),cds(2),'yo','markersize',20);
            hold off;
            axis off;
        end; % showRadialAverage
    drawnow;
    end; % for j=1:nIons
end; % for i-1:nAngs
%%
if showRadialAverage
%     figure(2*bigFigure);
%     plot(xVals,squeeze(lines(:,1,1,:)));
%     legend(legTxt);
%     title('Means');
    if bigFigure==2
        figText='Na_I_Map';
        maxY=1.5;
    else
        figText='Na_Map';
        maxY=1;
    end;
    figure(2*bigFigure+1);
    plot(xVals,squeeze(lines(:,2,1,:)),'linewidth',1);
    hold on;
    plot(xVals,0*xVals,'k-');
    hold off;
    axis([min(xVals) max(xVals) -.5 maxY]);
    
    
    legend(legTxt);
    title([figText ' Medians'],'interpreter','none');
    xlabel('Radius from atom center, Å');
    ylabel('Map density');
end;


return

% %%
%
% msk=fuzzymask(n,3,[45,20,45]);
%
%
%
% %
% fc=.5;
% mdiff=GaussFilt(m2cr,fc)-m1cr;
% figure(3);
% ShowSections(mdiff);
%
% %% Upsample the maps
%
%
%% locate the ions
% p=pdb2mat('SR_7_1_21_dimer_copy.pdb');
% [u,a,b]=unique(p1.element);
% ptrsNa=find(b==4);
% ptrsI=find(b==2);
% iPtrs=[ptrsNa(2:3); ptrsI(1)];
% iPtrs2=[ptrsNa(1); ptrsNa(4); ptrsI(2)];
%
% cts=n/2-1; % center shift: 288 -> 144, plus center is 73 not 72.
% iXs=p1.X(iPtrs)/s.pixA-cts;
% iYs=p1.Y(iPtrs)/s.pixA-cts;
% iZs=p1.Z(iPtrs)/s.pixA-cts;
% iXs2=p1.X(iPtrs2)/s.pixA-cts;
% iYs2=p1.Y(iPtrs2)/s.pixA-cts;
% iZs2=p1.Z(iPtrs2)/s.pixA-cts;
%
% %%
% % pdb coords have the center being
% figure(4);
% ct1=n/2+1;
% slice=sum(m1c(:,:,4:6+ct1),3);
% rSlice=rsRotateImage(slice,180);
% imags(slice+rSlice);
% hold on;
% plot(iXs,iYs,'yo');
% plot(iXs2,iYs2,'yo');
% % plot(73-iXs2+73,73-iYs2+73,'ro');
% plot(ct1,ct1,'y+');
% hold off;
%
% %% plot z densities
% zds=zeros(n,3);
% for i=1:3
%     zds(:,i)=m1c(round(iXs(i)),round(iYs(i)),:);
% end;
% figure(5);
% plot(zds(73:81,:));
%
%% Get the rotated ion coords

% upsample
us=4;
m1cr4=Downsample(m1cr,us*n);
ct1=us*n/2+1;
iXsr=us*p1r.X(iPtrs);  % zero-based
iYsr=us*p1r.Y(iPtrs);
iZsr=us*p1r.Z(iPtrs);
%%
for k=1:3
    %     sliceR=sum(m1cr4(:,:,round(iZsr(k))+ct1),3);
    sliceR=m1cr4(:,:,round(iZsr(k))+ct1);
    figure(3+k);
    imags(sliceR);
    hold on;
    plot(iXsr(k)+ct1,iYsr(k)+ct1,'yo');
    %     text(iXsr(i)+ct1+us,iYsr(i)+ct1+us,num2str(i),'color',[1 1 0]);
    text(iXsr(k)+ct1+us,iYsr(k)+ct1+us,num2str(k),'color',[1 1 0]);
    % plot(iXs2,iYs2,'yo');
    % plot(73-iXs2+73,73-iYs2+73,'ro');
    plot(ct1,ct1,'y+');
    hold off;
end;

%%
% plot z densities again
zR0=5;
zR=us*zR0; % we'll compute +- this z value.
zs=-zR0:1/us:zR0-1/us;
zds=zeros(2*zR,3);
for i=1:3
    zds(:,i)=m1cr4(round(iXsr(i)+ct1),round(iYsr(i)+ct1),round(iZsr(i)-zR+ct1):round(iZsr(i)+zR-1+ct1));
end;
figure(7);
plot(zs,zds);
legend(p1.element(iPtrs));

return
%%

m1crf=GaussFilt(m1cr,.2);
fsh=FourierShift(n*[1 1 1],[0 0 -5]);
m1crf=real(ifftn(fftn(m1crf).*fsh));

% Aligner
%  two variables: alpha, dz
iVals=[34 0 .2 .2 .5 .5];
iSteps=[.2 1 .05 .05 .1 .1];
iFlags=[0 0 1 1 1 1];

P1=Simplex('init',iVals,iSteps,iFlags);
iMax=80
for i=1:iMax
    P=P1;
    m1crf=rsRotateImage(m1c,P(1));
    fsh=FourierShift(n*[1 1 1],[0 0 P(2)]);
    m1crs=real(ifftn(fftn(m1crf).*fsh));
    m1crs2=P(5)*GaussFilt(m1crs,P(4))+P(6)*SharpFilt(m1crs,P(5));
    figure(1);
    ShowSections(m1crs2);
    figure(3);
    mdiff=(m1crs2-m2cr).*msk;
    mdiff=GaussHP(mdiff,.05);
    ShowSections(mdiff);
    title(i);
    err=mdiff(:)'*mdiff(:);
    disp([num2str([i P]) '   ' num2str(err)]);
    P1=Simplex(err);
    if i==iMax/2
        P1=Simplex('init',P,iSteps,iFlags);
    end;
    pause(0.1);
end;
