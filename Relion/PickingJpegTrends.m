% PickingJpegTrends.m
% Take the overlaps and intensity data and try to combine them.
cd('/Volumes/D257/Hideki/201228/05_1_1/Picking_8/');
load('allMis8_intens2.mat')
%  loads allMis, basePath
 
nm=numel(allMis);
 
oks=zeros(nm,20);
for i=1:nm
    oks(i,:)=allMis{i}.ok(1:20);
end;
 
plot(oks);

load frac

return

% %%
% % Get the fraction overlaps
% frac=zeros(nm,1);
% for i=1:nm
%     [m,frac(i)]=FractionOverlaps(allMis{i});
%     imags(m); title([i frac(i)]');
%     drawnow;
% end;
% % save frac frac
% return

%%
% look at trends between the intensity and frac

plot(oks(:,15),frac,'o');

oki=oks(:,15);
fri=frac;

figure(3);
indsOk=(oki>7.5 & oki<8.7) & fri<.25;
hist(oki(indsOk),7.5:.01:8.7);
figure(4);
hist(fri(indsOk),0:.001:.25);
figure(2);
plot(fri(indsOk),oki(indsOk),'.');
sum(indsOk)

hold on;
activeSet=oki>8.1 & oki<8.5 & fri<.06; % selects about 7500 micrographs
 plot(fri(activeSet),oki(activeSet),'k.');
 hold off;
 
 % assign frac values to mis
 
 for i=1:nm
     allMis{i}.ok(2)=frac(i);
     allMis{i}.ok(8)=activeSet(i);
     allMis{i}.active=activeSet(i);
 end;
 sum(activeSet)
 
 allMiName='allMis_intens+frac_7505';
 disp(['Saving ' allMiName]);
 
 save(allMiName, 'allMis');
 
