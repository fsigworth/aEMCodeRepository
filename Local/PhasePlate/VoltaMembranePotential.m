% VoltaMembranePotential.m
% Compute image effect of a membrane potential in a vesicle in the case of
% volta phase contrast

figure(1);
SetGrayscale;
jpegDir='/Users/fred/matlabWork/PhasePlate/Jpeg/';

% cutOn=300;  % A

pixA=4;
n=ceil(1000/pixA)*2;
nh=n/2-1;
xs=0:pixA:(n-1)*pixA;

% Create the membrane model
vLipid=1.6;  % inner potential in volts
V0=0.2;  % internal potential
sigmaE=.73; % phase shift (mrad) per V.A  200 kV
thk=50;
rise=8;
% Create the model, which is sampled in units of the original pixel size.
nm0=ceil((thk/2+rise*1.5)/pixA)*2+1;  % array for vesicle model; 60A nominal
membraneModel=fuzzymask(nm0,1,thk/pixA/2,rise/pixA)...
    *vLipid*pixA;  % units of V.A per voxel

as=10.^(2:.1:3);
as=[100 150 200 300];
na=numel(as);
cutOns=[100 200 300 500 1000 2000 2000];
cutOns=[400 800];
ncs=numel(cutOns);
iAs=zeros(na,ncs,2);
    
    for i=1:na % loop ovr radii
        a=as(i);
        
        nm1=ceil(a+thk)*2+1;
        sectModel=fuzzymask(nm1,1,a/pixA,rise/pixA)*pixA;
        interiorModel=sectModel((nm1+1)/2:nm1);
        nm2=numel(interiorModel);
        interior=VesicleFromModel(n,nm2/2,interiorModel)*sigmaE*V0;
        
        vesicle=VesicleFromModel(n,a/pixA,membraneModel)*sigmaE;
        
        
        %% filtering
        riseFraction=.25;
        linePlots=[];
for j=1:ncs
    cutOn=cutOns(j);
    
        
        fc=pixA/cutOn;
        filtVes=SharpHP(vesicle,fc,fc*riseFraction);
        filtInt=SharpHP(interior,fc,fc*riseFraction);
        
        subplot(221);
        imacs(xs,xs,filtVes);
        xlabel('X-coordinate, Å');
        title([num2str(cutOn) ' Å cut-on']);
        subplot(222);
        imacs(xs,xs,filtInt);
        title([num2str(a) 'Å radius, ' num2str(V0*1000) 'mV']);
        linePlots=[linePlots sect(filtVes) sect(filtInt)];
        subplot(223);
        plot([sect(filtVes) sect(filtInt)]);
        
        filtBoth=filtVes+filtInt;
        filtBoth(:,n/2+1:n)=filtVes(:,n/2+1:n);
        spNumber=10+mod(j,2)+4*round((j-1)/4)
        subplot(4,4,spNumber);
        imacs(Crop(filtBoth,n/2));
%         
%         
%         q1=filtVes(:)/1000;  % radians
%         q2=filtInt(:)/1000;
%         A=zeros(2,2);
%         A(1,1)=q1'*q1;
%         A(1,2)=q1'*q2;
%         A(2,1)=A(1,2);
%         A(2,2)=q2'*q2;
%         iA=inv(A);
%         % iA
%         E=[q1 q2]*iA;
%         e1=reshape(E(:,1),n,n);
%         e2=reshape(E(:,2),n,n);
%         subplot(224);
%         imacs(xs,xs,e2);
%         drawnow;
%         iAs(i,j,1)=iA(1,1);
%         iAs(i,j,2)=iA(2,2);
    end;
        subplot(223);
        plot(linePlots);
    pause
end;
%%
% subplot(223);
% loglog(as,sqrt(iAs(:,:,2)),'-');
% xlabel('Radius, Å');
% ylabel('Standard error');
% title('Em est error');
% legend(num2str(cutOns'));
% subplot(224);
% loglog(as,sqrt(iAs(:,:,1)),'-');
% xlabel('Radius, Å');
% ylabel('Standard error');
% title('Ves est error');
% legend(num2str(cutOns'));
% 
% name=sprintf('VesMbnPot%dCO.jpg',cutOn);
% set(gcf,'paperpositionmode','auto');
% print('-djpeg','-r150',[jpegDir name]);
% 
