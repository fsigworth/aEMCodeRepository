% ShowHemispheres.m


mapName='/Volumes/TetraData/Structures/AMPAR/3KG2map58.mrc';
[map pixA]=ReadEMFile(mapName);

angs=[0 80 0; 0 100 0
      0 60 0; 0 120 0
      0 40 0; 0 140 0
      0 20 0; 0 160 0];

% % Look at variations in gamma
% for i=1:8
%     angs(i,:)=[0 40 (i-1)*25];
% end;
% 

figure(1);
SetGrayscale;
projs=rsMakeTemplates(angs,map);
na=size(angs,1);
for i=1:na    
    subplot(floor((na+1)/2),2,i);
    imacs(projs(:,:,i));
    title(angs(i,2));
end;
drawnow;
