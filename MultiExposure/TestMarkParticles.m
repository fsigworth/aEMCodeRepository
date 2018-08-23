function mts=TestMarkParticles(mi,mts)
nim=size(mts,3);
for i=1:nim
    m=mts(:,:,i);
    mxm=max(m(:));
    np=size(mi.particle.picks,1);
    for j=1:np  % loop over particles
        flag=mi.particle.picks(j,3);
        if flag>=16 && flag <=33  % a real particle
            coords=mi.particle.picks(j,1:2)+1;
            vind=mi.particle.picks(j,4); % vesicle index
            if vind>0
                shCoords=coords+[mi.vesicle.shiftX(vind,i) mi.vesicle.shiftY(vind,i)];
            else
                shCoords=coords;
            end;
%             Draw a white box
            for ix=-1:1
                for iy=-1:1
            m(ix+round(shCoords(1)),iy+round(shCoords(2)))=mxm;  % mark a bright spot.
                end;
            end;
        end;
    end;
    mts(:,:,i)=m;
end;
