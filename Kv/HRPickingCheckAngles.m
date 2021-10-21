% HRPickingCheckAngles.m
% See if we are sampling angles with small enough steps

refAngles=[10 20 180]; % phi, theta, psi
shiftAngles=[0 0 0; 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
refDir='HRPicking/';
refName='tmMap.mrc';
zsh=[0 0 21]; % shift TM region to center

ds=3;
n=48
angleSteps=4; % Figure 3 A steps at a radius of 40 A
[map3,s]=ReadMRC([refDir refName]);
pixA=s.pixA*ds;
map=DownsampleGeneral(circshift(map3,zsh),n,1/ds);
ShowSections(map);
%%
aShift=0;
xyShift=0.33;
xyShifts=shiftAngles(:,1:2)*xyShift; % use the first two elements of shiftAngles
nsh=size(shiftAngles,1);
angs=aShift*shiftAngles+repmat(refAngles,nsh,1);
projs=rlMakeTemplates(angs,map,0,xyShifts);

ccs=zeros(n,n,nsh);
maxs=zeros(nsh,1);
fRef=conj(fftn(ifftshift(projs(:,:,1))));
figure(2);
for i=1:nsh
    ccs(:,:,i)=real(ifftn(fftn(projs(:,:,i)).*fRef));
    subplot(3,3,i);
    imags(ccs(:,:,i));
    maxs(i)=max2d(ccs(:,:,i));
    title(maxs(i)/maxs(1));
end;
drawnow;