sumImg=zeros(128,128);
nim=size(s,3);
for i=1:nim
    mi=si.mi{si.miIndex(i)};
    iso=mi.particle.picks(si.miParticle(i),7);
    if ~iso
        sumImg=sumImg+s(:,:,i);
    end;
end;
