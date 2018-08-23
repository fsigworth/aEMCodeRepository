function [mi,vesIndex,partIndex,vesStruct,picksArray]=rsVesicleIndFromSi(si,pIndex)
% Given the particle index, return information about its vesicle.
% mi is the relevant mi structure, vesIndex is the vesicle in that
% micrograph, partIndex is the particle index in that micrograph.
% vesStruct is simply mi.vesicle; picksArray is mi.particle.picks(partIndex,:).
% Thus you can get vesStruct.x(vesIndex), vesStruct.s(vesIndex,:) etc.
% and also picksArray=[ x y flag ves ... ] for the particle.

imi=si.miIndex(pIndex);
mi=si.mi{imi};
partIndex=si.miParticle(pIndex);
picksArray=mi.particle.picks(partIndex,:);
vesIndex=picksArray(4);
vesStruct=mi.vesicle;
