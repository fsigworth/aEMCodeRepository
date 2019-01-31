% rsPickerStats.m
% run after rsCountParticles.m

% MiLoadAll;
%%

infoDir='Info/';
    load([infoDir 'allMis.mat']);
    nEntries=numel(allMis);
nd=16;
minDef=1;
maxDef=2;
vars=[];
amps=[];
points=0;
ptrs=zeros(nEntries+1,1);
defOk=false(nEntries);

for i=1:nEntries
    mi=allMis{i};
    d=mi.ctf(1).defocus;
    if numel(mi.particle.picks)>0 && d>=minDef && d<maxDef
        defOk(i)=true;
        ptrs(i)=numel(amps)+1; % pointer into the huge amps array.
        amps=[amps; mi.particle.picks(:,5)];
        vars=[vars; mi.particle.picks(:,8)];
        points=points+1;
    end;
end;
ptrs(nEntries+1)=numel(amps)+1;
points

gAmpMin=3.7; % 2.9; %3;
gAmpMax=4.5; %3.8; %3.4;
gVarMin=0;
gVarMax=4;

nq=numel(amps);
oks=false(nq,1);
oks=amps>=gAmpMin & amps<gAmpMax & vars>=gVarMin & vars<gVarMax;
goods=sum(oks)



aMin=min(amps);
aMax=max(amps);
vMin=min(vars);
vMax=max(vars);
sclAmps=floor((nd-1)*(amps-aMin)/(aMax-aMin))+1;
sclVars=floor((nd-1)*(vars-vMin)/(vMax-vMin))+1;
img=zeros(nd,nd,'single');
for i=1:numel(sclAmps)
    img(sclVars(i),sclAmps(i))=img(sclVars(i),sclAmps(i))+1;
end;

figure(3);
subplot(221);
[h,b]=hist(amps,100);
bar(b,h);
xlabel('Amplitudes');
ylabel('Frequency');
title(['Var range ' num2str([gVarMin gVarMax],3)]);

subplot(222);
bar(b,cumsum(h,'reverse'));
% plot(vars,amps,'.');
xlabel('Amps');
ylabel('Cum histo');
% axis([vMin vMax aMin aMax]);
title(['Defocus ' num2str([minDef maxDef])]);

subplot(223);
hist(vars,100);
xlabel('Vars');
ylabel('Frequency');
title(['Amp range ' num2str([gAmpMin gAmpMax],3)]);

figure(3);
subplot(224);
k=.5;
k=.25;
fimg=GaussFilt(img,.1);
imac([vMin vMax],[aMin aMax],k*fimg);
axis xy
ylabel('Amp');
xlabel('Var');
colorbar;
title([num2str(points) ' images']);


return
% good region: var 5.5-8; amp 3-3.4;
%% Write modified info files.
allMis1=allMis;
for i=1:nEntries
    mi=allMis{i};
    mi1=mi;
    if defOk(i)
    partOks=oks(ptrs(i):ptrs(i+1)-1);
    mi1.particle.picks=mi.particle.picks(partOks,:);
    end;
    allMis1{i}=mi1;
end;
% disp('saving allMis1.mat');
% save([infoDir 'allMis1.mat'],'allMis1');
% 


