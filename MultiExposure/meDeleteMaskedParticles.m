% meDeleteMaskedParticles
% Find all the particles in the hole mask, and delete them

displayOn=0;
% load allMis_holes_i2_ov_cls.mat
nim=numel(allMis);
newAllMis=allMis;
totalDeleted=0;
for i=1:nim
    numDeleted=0;
    mi=allMis{i};
    if i==1
        maskSize=mi.mask(2).data.sz;
        M0=meGetImageScaling(mi.imageSize,maskSize);
        M1=inv(M0);
    end;
    smMask=meGetMask(mi,maskSize,2);
    if isfield(mi,'particle') && isfield(mi.particle,'picks') && size(mi.particle.picks,1)>0
        np=size(mi.particle.picks,1);
        dels=false(np,1);
        disMask=smMask;
        disMask(1,1)=true;
        disMask(2,1)=false;
        if displayOn
            imags(disMask);
            hold on;
        end;
        for j=1:np
            xy=round(M1*[mi.particle.picks(j,1:2) 1]');
            if displayOn
                plot(xy(1),xy(2),'y+');
            end;
            xz=max([1 1],min(maskSize,xy(1:2)'));
            dels(j)= ~smMask(xz(1),xz(2));
        end;
        hold off;
        drawnow;
        numDeleted=sum(dels);
        if numDeleted>0
            mi.particle.picks=mi.particle.picks(~dels,:);
            newAllMis{i}=mi;
            totalDeleted=totalDeleted+numDeleted;
        end;
    end;
    fprintf('%5d %4d %4d %s\n',i,numDeleted,np, mi.baseFilename);
end;
totalDeleted
% have to copy newAllMis to allMis and save it, manually.
