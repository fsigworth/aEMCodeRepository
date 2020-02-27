function [nParts,partsPerVes,partsPerArea]=rspGetParticleStats(mi);
% Compute the number of particles, particles per total vesicles, and
%  particles per total r^2 of all vesicles up to the largest vesicle on
%  which particles were found.

%  Default return values
nParts=0;
partsPerVes=0;
partsPerArea=0;

if numel(mi.particle.picks)>0
    flags=mi.particle.picks(:,3);
    parts=(flags>=16 & flags<48);
    if sum(parts)>0 % we have some particles to analyze
        vesInds=mi.particle.picks(parts,4);
        radii=real(mi.vesicle.r(vesInds,1));
        maxRadius=max(radii);
        allRadii=real(mi.vesicle.r(:,1));
        goodVes=mi.vesicle.ok(:,2);
        okVes=goodVes & allRadii<=maxRadius;
        nParts=sum(parts);
        partsPerVes=nParts/sum(goodVes);
        partsPerArea=1e6*nParts/sum(allRadii(okVes).^2);
    end;
end;
