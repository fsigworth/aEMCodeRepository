% IonPeaksFromMap.m
%  Read Silvia's map and examine the vicinity of the ion densities.

% cd('/Users/fred/EMWork/Silvia');
% [map,s]=ReadMRC('relion678_150.mrc');
% 
% d=dir;
% [coords,rots,type]=ReadPDB(d(3).name);
% iCoords=coords(:,end-5:end); % The last 6 coords are ions.
iNames={'Na1C' 'Na2C' 'Na3C' 'Na1D' 'Na4D' 'I1E'};

for i=1:6
cp=iCoords(:,i)/s.pixA+1;
ms=ExtractVolume(map,round(cp),32);
ShowSections(ms);
subplot(339);
r=Radial3(map,cp);
plot(r(1:32),'.-','linewidth',1);
xlabel('Radial distance, A');
ylabel('Shell average map value');
mysubplot(331);
title(iNames{i});
drawnow;
print('-djpeg',[iNames{i} '.jpg']);

% % make a fake atom
% m6=fuzzymask(256,3,3,1,cp);

end;

%%
% Compute the vector between the Na2C and I1E coordinates, and try to align
% the peaks.  note that we still need to make good fractional shifts...

vecNI=iCoords(:,2)-iCoords(:,6);
% we find that EulerMatrix(.845,.12,0)*vecNI is a vector along X direction,
% approximately [-5 0 0].
angs=[.845 .12 0];
invAngs=[0 -.12 -.845]; % ERotate by this to make the vector in X.

% Get the vicinity of I1E
cp=iCoords(:,6)/s.pixA+1;
ms=ExtractVolume(map,round(cp),64);
ctr=33*[1 1 1];
figure(1);
ShowSections(Crop(ms,32));
r=Radial3(ms,ctr);
plot(r,'.-','linewidth',1);
xlabel('Radial distance, A');
ylabel('Shell average map value');


msr=ERotate3(ms,angs);
figure(2);
ShowSections(Crop(msr,32));
r=Radial3(msr,ctr);
plot(r,'.-','linewidth',1);
xlabel('Radial distance, A');
ylabel('Shell average map value');


%% put the map into standard position
 angs=[35+90 90 90]*pi/180;
% mapd=Downsample(Crop(map,128),256);
% mapdr=ERotate3(mapd,angs);
% ShowSections(mapdr);
clf;
rCoords=EulerMatrix(angs)*coords/1.068;
plot3(rCoords(1,:),rCoords(2,:),rCoords(3,:),'.');




