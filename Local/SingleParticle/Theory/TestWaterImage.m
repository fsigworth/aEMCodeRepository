% TestWaterImage

cd /Users/fred/EMWork/waterbox/pdbs

disp('Reading pdb file');
[coords, types]=ReadPDBAtoms2('wb99.pdb');
disp('done.');
%%
nAtoms=size(coords,2);
oxy=false(nAtoms,1);
for i=1:nAtoms
    oxy(i)=strncmp(types(i,:),'OH2',3);
end;
oxyPos=coords(:,oxy)';
nOxy=size(oxyPos,1);
mins=min(oxyPos);
pos=oxyPos-repmat(mins,nOxy,1);
maxs=max(pos);
n=NextNiceNumber(maxs(1:2));
os=4;
imgOs=zeros(n*os,'single');

for i=1:nOxy
    ix=floor(os*pos(i,1))+1;
    iy=floor(os*pos(i,2))+1;
    imgOs(ix,iy)=imgOs(ix,iy)+1;
end;
imgOs=Crop(imgOs,768);
img=Downsample(imgOs,192);

figure(1);
SetGrayscale;
imacs(img);
