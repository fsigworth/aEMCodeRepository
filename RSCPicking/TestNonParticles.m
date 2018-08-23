% TestNonParticles
nt=40;
r0=16;
npAngleList=rsListHemisphereAngles(32, 12);
na=size(npAngleList,1)
nr=0;
    rds=[5 8 16 32];
for j=1:numel(rds);
    rd=rds(j);
r1=rd+r0;
subplot(223);
% small disc
for i=1:na
    r2=r1*npAngleList(i,2)/90;  % don't take sine.
    x=nt/2+1+r2*sind(npAngleList(i,1));
    y=nt/2+1+r2*cosd(npAngleList(i,1));
    im=fuzzymask(nt,2,rd,2,[x y]).*fuzzymask(nt,2,r0,3);
    nr=nr+1;
    npRefs(:,:,i,j)=im;
%     imacs(im);
%     drawnow;
end;
end;
nrds=numel(rds);