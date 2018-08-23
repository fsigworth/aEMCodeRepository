% ScanCheckCTFs.m

field='Defocus, µm';
nv=1;  % pick up this many values

% Have the user select some mi files
[fname, pa]=uigetfile('*mi.*','Select mi files to examine','multiselect','on');
if ~iscell(fname)
    fname={fname};
end;
[rootPath,infoPath]=ParsePath(pa);
%%

cd(pa);
nim=numel(fname);
val=zeros(nim,nv);
for i=1:nim
    %%
    mi=ReadMiFile(fname{i});
    mi.basePath=rootPath;
            val(i,1)=mi.ctf(1).defocus;
% Fit the CTF
m=meReadMergedImage(mi);
fitPars=mi.ctf(1);
fitPars.defocus=1:.1:7;
fitPars.deltadef=-.2:.05:.2;
fitPars.theta=0:pi/10:pi;
P1=FitCTF(m,fitPars,mi.pixA*2,8,15);
P2=FitCTF(m,fitPars,mi.pixA*2,10,30);

% opts=struct;
% opts.maxRes=8;
% opts.minRes=20;
% opts.kV=300;
% P3=meFitCTF(m,mi.pixA*2,4,0,opts);
% P3
title(fname{i},'interpreter','none');
val(i,2)=P1.defocus;
val(i,3)=P2.defocus;
if any(abs(val(i,2:3)-val(i,1))>.1)
    str='****';
else
    str='';
end;
disp([fname{i} '  ' num2str(val(i,:)) '  ' str]);
end;
return


%% Code to re-do CTF fitting and MergeImages
cd /Users/fred/EMWork/Hideki/151216
load CTFErrs.mat  %% get val, fname
%%
active=abs(val(:,1)-val(:,3))>.1;
numToProcess=sum(active)
namesToProcess=fname(active);
fname0=fname;
infoNames=namesToProcess;
for i=1:numToProcess
    infoNames{i}=['Info/' namesToProcess{i}];
end;

mpars.overwrite=1;
ctopts.defoci=1:.2:10;
ctopts.maxRes=8;
ctopts.minRes=30;
mpars.ctfOptions=ctopts;

% MergeImages(infoNames,mpars);
%%
fname=namesToProcess;
infoPath='Info/';
rootPath=AddSlash(pwd);
rsRefineVesicleFits;






