% rtDefocusSeriesEval
% From all the .mrc files in a directory, accumulate the Gctf results.

% cd '/Volumes/Drobo4/201909X_RIKEN/190902/RELION/CtfFind/job018/Movies'
cd '/Volumes/Drobo4/201909X_RIKEN/190902/RELION/CtfFind/job017/Split1536' % 3x3
% cd '/Volumes/Drobo4/201909X_RIKEN/190902/RELION/CtfFind/job023/Movies'
startIndex=4;
multiRun=5;
if multiRun
    startIndex=10;
end;

% cd '/Volumes/Drobo4/201909X_RIKEN/190902/RELION/CtfFind/job017/Split1536'
% startIndex=1;
% multiRun=5;

% cd '~/EMWork/Hideki/201909XX/190902/Relion/CtfFind/job009/Movies'
% cd '~/EMWork/Hideki/201909XX/190902/Relion/CtfFind/job017/Split1536'

writeJpeg=0;

d=dir;
ext='.mrc';
names={};
ni=0;
for i=1:numel(d)
    if ~d(i).isdir && strndcmp(d(i).name,ext)
        ni=ni+1;
        names{ni}=d(i).name;
    end;
end;
names(1:startIndex-1)=[];
ni=ni-startIndex+1;

for i=1:ni
%     names{i}
  [epaVal,ctfVal,ctfImg]=rtReadGctfLogs('',names{i});
  if i==1
      epaVals=epaVal;
       ctfVals=ctfVal;
       ctfImgs=ctfImg;
  else
      epaVals(i)=epaVal;
      ctfVals(i)=ctfVal;
      ctfImgs(:,:,i)=ctfImg;
  end;
end;
%%
ndp=floor(size(epaVal.epaRaw,1)/2);
dp=1:ndp;

freqs=1./epaVals(1).resolution(dp);
legends={};
% nepa=numel(epaVals(1).ctfSim);
nepa=ndp;


for mOffset=1:8
if multiRun
iOffset=mOffset*9;
ni=9;
else
    iOffset=0;
    mOffset=0;
end;

epaRaw=zeros(ndp,ni);
epaBkgSub=zeros(ndp,ni);
epaCC=zeros(ndp,ni);
epaSim=zeros(ndp,ni);
cc=zeros(ndp,ni);
defocus=zeros(ni,1);
astig=zeros(ni,1);
legends=cell(ni,1);

for i=1:ni
    j=i+iOffset;
%     if i==1
%         epaRaw=epaVals(1).epaRaw(dp);
%         epaBkgSub=epaVals(1).epaBkgSub(dp);
%         epaCC=epaVals(1).ccc(dp);
%         defocus=(ctfVals(1).Defocus_U+ctfVals(1).Defocus_V)/2e4;
%         astig=(ctfVals(i).Defocus_U-ctfVals(i).Defocus_V)/2e4;
%         cc=(epaVals(1).ctfSim(dp).^2-.5).*epaVals(1).epaBkgSub(dp);
% %         zeros=fin(epaVals(i).epaSim)
%     else
        epaRaw(:,i)=epaVals(j).epaRaw(dp);
        epaBkgSub(:,i)=epaVals(j).epaBkgSub(dp);
        epaCC(:,i)=epaVals(j).ccc(dp);
        epaSim(:,i)=epaVals(j).ctfSim(dp);
        cc(:,i)=(epaVals(j).ctfSim(dp).^2-.5).*epaVals(j).epaBkgSub(dp);
        defocus(i)=(ctfVals(j).Defocus_U+ctfVals(j).Defocus_V)/2e4;
        astig(i)=(ctfVals(j).Defocus_U-ctfVals(j).Defocus_V)/2e4;
%     end;
%     legends{i}=num2str([defocus(i) astig(i)],3);
    legends{i}=[num2str(defocus(i),4) '   ' num2str(astig(i),2) '    ' ...
        num2str(ctfVals(j).Angle,3)];
    
end;
% Compute equiphase models

cPars=struct;
cPars.lambda=EWavelength(200);
cPars.Cs=2.7;
cPars.B=0;
cPars.alpha=.1;
cPars.pixA=.67;
nf=numel(freqs);
model1D=zeros(nf,ni);
chi1D=zeros(nf,ni);
for i=1:ni
    j=i+iOffset;
    defU=ctfVals(j).Defocus_U;
    defV=ctfVals(j).Defocus_V;
    angle=ctfVals(j).Angle;
    cPars.defocus=1e-4*(defU*cosd(angle)^2+defV*sind(angle)^2);
    [model1D(:,i),chi1D(:,i)]=ContrastTransfer(freqs,cPars);
end;


%Compute correlation coefficients
% Do this simultaneously for all the traces
m1=model1D.^2-.5;
m2=m1.^2;



cycles=round(-chi1D);

% average baseline over cycles
epaBase=zeros(nf,ni);
bsEpa=zeros(nf,ni);
for i=1:ni
    [c,ia,ic]=unique(cycles(:,i));
    [sumE1,norm]=WtHist(ic,epaBkgSub,max(ic));
    epaBase(:,i)=sumE1(ic)./norm(ic);
    bsEpa(:,i)=epaBkgSub(:,i)-0*epaBase(:,i);
end;


% e2=epaBkgSub.^2;
% cc=m1.*epaBkgSub;
e2=bsEpa.^2;
cc=m1.*bsEpa;

% average over each cycle
exCC=zeros(nf,ni);
for i=1:ni
    [c,ia,ic]=unique(cycles(:,i));
    [sumCC,norm]=WtHist(ic,cc(:,i),max(ic));
    sumM2=WtHist(ic,m2,max(ic));
    sumE2=WtHist(ic,e2,max(ic));
    avgCC=min(1,sumCC./(sqrt(sumM2.*sumE2)));
    exCC(:,i)=avgCC(ic);
    
end;

% figure(2);
%     plot(exCC);
% 
lw=2;

shift=0:-1:1-ni; % shift for display

% ------main plot-------
figure(1);
mShift=repmat(shift,dp(end),1);
% plot(freqs(dp),epaRaw(dp,1:end)+mShift);


set(gca,'ColorOrderIndex',1);
% plot(freqs(dp),(model1D(dp,:).^2-.5)*0.3+mShift);
plot(freqs(dp),(epaSim(dp,:).^2)*0.1+mShift,'-','color',[.7 .7 .7] ,'linewidth',lw);
hold on;


set(gca,'ColorOrderIndex',1);
% plot(freqs(dp),epaBkgSub(dp,1:end)+mShift);
plot(freqs(dp),2*bsEpa(dp,1:end)+mShift,'linewidth',lw); % epaBkgSub


set(gca,'ColorOrderIndex',1);
% plot(freqs(dp),epaBase(dp,:)+mShift,'k-');
% plot(freqs(dp),cc(dp,:)+mShift,'k-');
plot(freqs(dp),epaCC(dp,:)/2+mShift,'k-');

set(gca,'ColorOrderIndex',1);
% plot(freqs(dp),exCC(dp,:)*0.5+mShift,'-');

plot(freqs(dp),mShift,'-','color',[.5 .5 .5]);
plot(freqs(dp),mShift+.5,'--','color',[.5 .5 .5]);

for i=1:ni
% text(.01,shift(i),legends{i},'fontsize',16);
str=[names{i} '  ' legends{i}];
text(freqs(round(ndp*.6)),shift(i)+.5,str,'fontsize',14,'interpreter','none',...
    'verticalalignment','bottom');
end
hold off;
% legend(legends(1:end));

if writeJpeg
    
    if multiRun
        str=num2str(mOffset);
    else
        str='All';
    end;
    
    print(['Fig' str '.jpg'],'-djpeg');
end;
if multiRun
    pause(2);
end;
end;


return
%%



figure(3);
for i=1:ni
    
imags(ctfImgs(:,:,i));
title(i);
pause(0.5)
end;

%% compute 1D spectra from the micrographs
sp1s=0;
shft=.1;
for i=1:ni
    disp(names{i});
    [m,s]=ReadMRC(names{i});
    sp1=RadialPowerSpectrum(m);

    if i==1
        sp1s=sp1;
    else
        sp1s(:,i)=sp1*shft;
    shft=shft/10;
    end;
end;
%%
nfp=size(sp1s,1);
dp=1:nfp/2;
f1s=dp/(2*nfp*s.pixA);
figure(3);
semilogy(f1s,sp1s(dp,:));
axis([0 inf 1e-4 1e6]);
legend(legends,'fontsize',18);  
    
    
    

