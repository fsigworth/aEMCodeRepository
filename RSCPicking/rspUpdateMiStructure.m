function mi=rspUpdateMiStructure(mi)
if ~isfield(mi.particle,'picks') % old-fashioned particle data
    np=numel(mi.particle.x);
    coords=single(zeros(np,8));
    for i=1:np
        coords(i,:)=[mi.particle.x(i) mi.particle.y(i) ...
            mi.particle.type(i) mi.particle.vesicle(i) 0 0];
    end;
    mi.particle.picks=coords;
end;
if isfield(mi.particle,'x')
    mi.particle=rmfield(mi.particle,{'x' 'y' 'vesicle' 'type' 'quality'});
end;

