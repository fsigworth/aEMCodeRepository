% rlSortIOParticles.m
% Given a particles.star file containing the .vesiclePsi field,
% compare that with rlnAnglePsi to determine whether a given particle is
% inside-out or right side out. Plot a histogram of the angle differences.
% Then create two new .star files, xxxx_iso.star and XXXX_rso.star where
% the particles have been sorted.

doWrite=1;
minTilt=20; % skip particles where tilt <20 or >160
inStarName='RSC/Class3Dj055_data_psi.star'
outStarBasename='RSC/Class3Dj055_data'
suffixes={'rso.star' 'iso.star'};

[nms,dats]=ReadStarFile(inStarName);

pts=dats{2}; % particles struct
tilts=pts.rlnAngleTilt;
okTilts=(tilts>minTilt & tilts<180-minTilt);
psis=pts.rlnAnglePsi;
psis=mod(psis,360);
vesPsi=-pts.vesiclePsi;
vesPsi=mod(vesPsi,360);

subplot(1,2,1);
plot(vesPsi,psis,'.','markersize',1);
xlabel('Vesicle psi');
ylabel('Particle psi');
subplot(1,2,2);
psiDiff=mod(vesPsi-psis+90,360);
hist(psiDiff,100);

rso=(psiDiff>180) & okTilts;
iso=(psiDiff<=180) & okTilts;
totalRSO=sum(rso);
totalISO=sum(iso);

txt=['No. of particles: ISO ' num2str(totalISO) '  RSO ' num2str(totalRSO)];
title(txt)
disp(txt);
disp([num2str(sum(~okTilts)) ' particle tilts rejected.']);

if doWrite
    
    flags=rso;
    dirTxt='rso';
    for i=1:2
        outName=[outStarBasename '_' dirTxt '.star'];
        disp(['Writing ' outName]);
        hdr=['# ver 30001 ' dirTxt ' particles from ' inStarName ' by rlSortIOParticles.m'];
        %     disp(hdr);
        WriteStarFile(nms,dats,outName,hdr,{[] flags});
        flags=iso;
        dirTxt='iso';
    end;
end;

disp(' done.');

