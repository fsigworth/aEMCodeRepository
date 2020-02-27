% GratingWav
% Animation of undiffracted waves

makeVideo=1;

b1=60; % left/lower border
b2=30; % other border
bs=b1+b2;
nx=1920-bs;
nz=1080-bs;
nt=100;

lambda=.02;
dz=lambda/10.5;
omega=.15;

rSeed=rng('default');

vName='UndiffWaves';

figure(2);
clf;
h=gcf;
h.Color=[.2 .2 .2];

pos=h.Position;
h.Position=[pos(1:2) nx+bs nz+bs];

if makeVideo
    v=VideoWriter([vName '.avi']);
            v.FrameRate=30;
    open(v);
    disp(['Making movie file: ' vName '.avi']);
end;

ha=axes;
ha.Units='pixels';
ha.Position=[b1 b1 nx nz];

zVals=(1:nz)*dz;
pars=struct;
pars.x=[1 nx]*dz;
pars.y=[1 nz]*dz;
pars.scl=1;
pars.sat=1.5;

rx=nx;
rz=nz;
fuzz=1;
envMag=.5;
envVel=40; % 40/1100 frames to move 2A means .06A/frame.
            % 3e8 m/s -> 3e18 A/s. 1 frame = .06/3e18 = 2e-20 s
rng(rSeed);
tmax=300;
t0=[0 60 180 240]';
xs=[0 nx nx*1.5 nx/2]';
orgs=[xs -2*rz-envVel*t0];
% orgs=[0 2*rz; nx 2*rz-envVel*t0(2); 
%     nx*1.5 2*rz-envVel*t0(3); nx/2 2*rz-envVel*t0(4)];
for t=0:10:tmax
                env=zeros(nx,nz,'single');
        for i=1:size(orgs,1)
%             if t>=t0(i)
                offset=orgs(i,:)+[0 1]*envVel*t;
                env=env+fuzzymask([nx nz],2,[rx rz],[rx rz]*fuzz, ...
                        offset);
%             end;
        end;
%         env=min(1,env);
        
        psi1=exp(1i*2*pi*(zVals/lambda-omega*t)); % undiffracted wave
psi=repmat(psi1,nx,1).*env;

imacx2(psi,1,pars);
axis ij
ha.FontSize=18;
ha.XColor=[1 1 1];
ha.YColor=[1 1 1];

xlabel('X, angstroms');
ylabel('Z, angstroms');
title([num2str(.01*t,'%6.2f'),' ns'],'FontSize',18,'Color',[1 1 1]);
drawnow;
    if makeVideo
    f=getframe(gcf);
    writeVideo(v,f);
    end;
end;

if makeVideo
    close(v);
end;
title('');
print('UndiffWave.jpg');

return

%% Figure with object
% after running GratingFresnel3 with lambda=.1, zBlock=.2;
% there dz=lambda/12.  So we plot more (nzn) points.
figure(2);
b1=70; % left/lower border
b2=20; % other border
bs=b1+b2;
ha.Position=[b1 b1 nx nz];
nx=1920-bs;
nz=1080-bs;
nzn=nz;
frac=0.7;
nz1=frac*nzn;
nz2=(1-frac)*nzn;
pars.x=xCal*10;
yExtent=nzn*.02/10;
pars.y=[-frac 1-frac]*yExtent;
psid1=psid;
psid1=psid1-.8*real(psid1);
imacx2(psid1(:,round(zOrigi-nz1):round(zOrigi+nz2)),1,pars);
% axis ij;
ha.FontSize=18;
ha.XColor=[1 1 1];
ha.YColor=[1 1 1];

xlabel('X, angstroms');
ylabel('Z, angstroms');

print('GratingWavesObject.jpg','-djpeg');





