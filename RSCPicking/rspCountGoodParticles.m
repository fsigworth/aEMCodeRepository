function count=rspCountGoodParticles(picks,dis)
minFlag=16;
maxFlag=47;
p=picks(:,:,3);
goods=p>=minFlag & p<=maxFlag;
count=[1 1]*sum(goods(:));
if dis.classParticlesMode && numel(dis.goodClasses)>0
    cls=dis.classes(:,:,1);
    flags=any(cls(:)==dis.goodClasses,2);
    count(1)=sum(goods(:).*flags);
% count=sum(flags(:)>=minFlag & flags(:) <= maxFlag);
end;


