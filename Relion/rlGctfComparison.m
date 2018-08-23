% rlGctfComparison.m

[bNames,q,ok]=ReadStarFile('Gctf2/micrographs_all_gctf.star');
q=q{1};

qDefs=(q.rlnDefocusU+q.rlnDefocusV)/2e4;
qDeltaDefs=(q.rlnDefocusU-q.rlnDefocusV)/2e4;
qThetas=q.rlnDefocusAngle;

load('Info/allMis.mat');
nmi=numel(allMis);

mDefs=zeros(nmi,1);
mDeltaDefs=zeros(nmi,1);
mThetas=zeros(nmi,1);
for i=1:nmi
    ctf=allMis{i}.ctf(1);
    mDefs(i)=ctf.defocus;
    mDeltaDefs(i)=ctf.deltadef;
    mThetas(i)=ctf.theta*180/pi;
end;

plot([mDefs qDefs]);

