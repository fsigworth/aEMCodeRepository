% MakeATPSynthaseMap.m

cd('/Users/fred/Structures/ATP Synthase')
[coords,types]=ReadPDBAtoms('5fik.pdb');
[composite,protein]=SolventAndProteinDensity(coords,types);

%%
% ShowSections(composite);
% ShowSections(protein);
solvent=composite-protein;
s0=solvent(1);
protMask=1-solvent/s0;  % 1 for protein region, 0 for solvent.
ShowSections(protMask);

%% Density modification in membrane region
n=size(composite,1)
ctr=ceil((n+1)/2)
mbnZ=80;
mbnCtr=48;
mbnAmp=2;
pixA=4;
% ShowSections(composite-protein,[ctr ctr mbnZ]);

subMap=composite-s0;
subMap(:,:,1:mbnZ)=(composite(:,:,1:mbnZ)-s0)-mbnAmp*protMask(:,:,1:mbnZ);

% ShowSections(GaussFilt(subMap,.1),[ctr ctr mbnCtr]);

mbnOffsetA=mbnCtr-ctr  % a negative number

% Downsampling

map=DownsampleGeneral(subMap,NextNiceNumber(n/pixA),1/pixA)*pixA;

ShowSections(map);
whos map
WriteMRC(map,pixA,'ATPSynthSub4A.mrc');
