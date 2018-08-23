% rlAddMembraneFromModel.m


% Put up a file selector for *si.mat,
disp('Getting an si file.');
[siName, siPath]=uigetfile('*si.mat','Select si file');
if isnumeric(siPath)  % user hasn't clicked Cancel
    return;
end;
cd(siPath);
disp(['Loading ' siName]);
si=load(siName);
si=si.si;
nim=numel(si.miIndex);
miIndex=si.miIndex(1);  % find the mi of the first particle
mi=si.mi{miIndex};

msk=ReadMRC('KvMask.mrc');

disp('Getting a subtracted map');

[vName, vPath]=uigetfile('*.mrc','Select map file');
if isnumeric(siPath)  % user hasn't clicked Cancel
    return;
end;
disp(['loading ' vPath vName]);
[vol,s]=ReadMRC([vPath vName]);
%%
ds=s.pixA/mi.pixA;
n=size(vol,1);

vds=meDownsampleVesicleModel(mi.vesicleModel,ds);
nds=numel(vds);
vCtr=(nds-1)/2;
ctr=ceil((n+1)/2);
r0=0.3*n;

vShift=round(si.mbnOffset);
vMembrane=zeros(n,n,n,'single');
for iz=1:nds
    vMembrane(:,:,vShift+ctr+iz-vCtr-1)=fuzzymask(n,2,r0,r0/10)*vds(iz);
end;

mScale=.03;

volM=vol+mScale*vMembrane;
ShowSections(volM,[ctr+4 ctr+4 ctr],45);

[pa,nm,ex]=fileparts([vPath vName]);
vNameRestored=[AddSlash(pa) nm 'MAdd' ex]
WriteMRC(volM,s.pixA,vNameRestored);

% 
% plot(vds);

