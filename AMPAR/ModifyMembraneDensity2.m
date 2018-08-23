% ModifyMembraneDensity2.m
% Reduce the TM density of the Kv channel to half.
pdbName='~/Structures/kv1.2/Kv1.2Tetramer.pdb';
miName='~/Structures/kv1.2/KvMembraneModel_mi.txt';
[coords, types]=ReadPDBAtoms(pdbName);
mi=ReadMiFile(miName);
tic
[totalDens, protDens]=SolventAndProteinDensity(coords, types);
toc
solDens=totalDens-protDens-totalDens(1); % subtract the background water dens.
%  protein will be at -4.95V.

%%
n=size(totalDens,1);  % it's 183
n1=96;

mctr=132; % membrane center is at y=132.
md0=mi.vesicleModel;
% pixA is 1.247
md=real(ifft(Cropo(fft(md0),71)));
modDens=totalDens-totalDens(1);
for iy=1:71
    y=mctr+iy-36;
    modDens(:,y,:)=modDens(:,y,:)+solDens(:,y,:)*md(iy)*0.33;
end;
% modDens=shiftdim(modDens,2);

modDens2=GaussFilt(Downsample(Crop(modDens,216),108),.2)*2;
modDens2=circshift(modDens2,[0 4 0]);
figure(4);
% ShowSections2(GaussFilt(modDens2-sim,.05));
q=shiftdim(GaussFilt(modDens2,.05)+GaussFilt(4*randn(108,108,108),.1),2);
q=GaussFilt(shiftdim(modDens2,2),.1);
ShowSections2(Crop(q,112),[],45);

sim=modDens2;
save('~/Structures/kv1.2/KvMapMbnSub.mat','sim');


return
%%
% Compare to our previous reconstructed image

% modDens1=DownsampleGeneral(modDens,96,1/(2*mi.pixA));

s=load('/Users/fred/Structures/kv1.2/KvMap.mat');
figure(1);
ShowSections2(shiftdim(Crop(s.sim,112),2),[],45);

% 
% figure(1);
% ShowSections2(GaussFilt(SharpFilt(modDens1,.12)+1.5*randn(n1,n1,n1),.17));

return
%%
map1Name='~/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon96j/mrc/i19av01.mrc';
map2Name='~/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon96j/mrc/i19bv01.mrc';
expMap=ReadMRC(map1Name)+ReadMRC(map2Name);
figure(2);
ShowSections2(expMap,[],45);



%%
% We use solDens as a mask.  

n1=size(TotalDens,1);
n0=2*ceil(n1/(2*pixA));  % size of final map, even to avoid possibe bugs

pixA=2.9;  % angstroms per pixel
npix=64;


cd('/Users/fred/aEMCodeRepository/AMPAR')
[m,pixA]=ReadEMFile('KvMap.mrc');

m=shiftdim(m,2);

m(:,:,65:100)=m(:,:,65:100)*.5;
% m(:,:,1:45)=0;
ShowSections2(m);

mbnOffsetA=48;
pixA=2;

map=m;

save KvMap.mat map pixA mbnOffsetA
