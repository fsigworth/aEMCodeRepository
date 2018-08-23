% MakeCTFGroups.m
% Open an si.mat file and cluster the CTF functions.  Create new fields in
% the si structure and store it.  The new fields:
% si.ctfGroupIndex : nim x 1    % index of group for each particle
% si.ctfGroups: n x n x ngroups % average CTF for each group.
nclasses=5;  % no. of defocus groups to make
nfactors=20;
freqExponent=.5;  % freq weighting of ctfs before clustering.  this is
% intended to line up the higher-frequency ripples.
% 
% [fname, pa]=uigetfile('*si.mat','Select a stack info file *si.mat');
% if isnumeric(pa) % File selection cancelled
%     return
% end;
% cd(pa);
% 
% load(fname);

%%
% Find the CTFs that aren't used
nmi=numel(si.mi);
h=hist(single(si.miIndex),1:nmi);
% There are some micrographs (and ctfs) that aren't used for any particles.
nullInds=find(h==0);
% Create the forward mapping indMap
% originalCtfIndex = indMap( activeCtfIndex)
indMap=(1:nmi)';
indMap(nullInds)=[];
% Create the reverse mapping from original CTFs to active CTFs
% activeCtfIndex = revMap(originalCtfIndex)
revMap=(1:nmi)';
for j=nullInds
    revMap(j+1:end)=revMap(j+1:end)-1;
end;
activeCtfs=si.ctfs(:,:,indMap);
% ctfs(:,:,nullInds)=[];
n=size(activeCtfs,1);
nctfs=size(activeCtfs,3);

%% Create defocus groups by clustering

boost0=Radius(n).^freqExponent;
boost=repmat(boost0,[1 1 nctfs]);

[means0, nm, inds, fvecs]=Classifier(activeCtfs.*boost,nclasses,nfactors,1);  % 20 factors, no decimation
%
means0=means0./repmat(boost0,[1 1 nclasses]);
figure(3);
ImagicDisplay1(activeCtfs);
figure(4);
ImagicDisplay1(means0);

nim=numel(si.miIndex);
si.ctfGroupIndex=zeros(nim,1,'int16');
si.ctfGroupIndex=uint16(inds(revMap(si.miIndex)));
si.ctfGroups=means0;

%%
% Compare the 1D ctfs with the class means
figure(1);
nr=floor(sqrt(nclasses));
nc=ceil(nclasses/nr);
for i=1:nclasses
    subplot(nr,nc,i);
    plot(sectr(means0(:,:,i)),'linewidth',2);
    hold on
    plot(sectr(activeCtfs(:,:,inds==i)));
    hold off;
end;

save(fname,'si');  % store the modified structure.
disp(['Saved the modified file ' fname]);
