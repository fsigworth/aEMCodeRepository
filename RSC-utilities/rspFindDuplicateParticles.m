function [inds,flags]=rspFindDuplicateParticles(mi)
% Search the particl.picks fields of an mi structure for particle locations
% closer than epsi.  If found, return the offending particle tuples and
% their flag (picks(:,3)) values in cell arrays.

epsi=4;
skipNullFlags=1;

inds={};
flags={};
nSets=0;
if isfield(mi,'particle') && isfield(mi.particle,'picks') && size(mi.particle.picks,1)>1
    np=size(mi.particle.picks,1);
    xs=mi.particle.picks(:,1);
    ys=mi.particle.picks(:,2);
    fs=mi.particle.picks(:,3);
    blanks=zeros(np,1);
    if skipNullFlags
        blanks(fs<16)=inf;  % not particles
        blanks(fs>=48)=inf; % not particles
    end;
    for i=1:np
        blanks(i)=inf;
        dist=sqrt((xs-xs(i)).^2+(ys-ys(i)).^2);
        dist(i)=inf;
        ind=find((dist+blanks)<=epsi);
        if numel(ind)>0
            blanks(ind)=inf;
            nSets=nSets+1;
            inds{nSets,1}=[i ind];
            flags{nSets,1}=fs(inds{nSets,1})';
        end;
    end;
end;
