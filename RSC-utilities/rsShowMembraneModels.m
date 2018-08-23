% ShowMembraneModels

    [si, imgs, names]=reLoadStackFiles;
    disp(names{end});
    
    nmi=numel(si.mi);
    model=zeros(0,nmi);
    j=0;
    for i=1:7:nmi
        if numel(si.mi{i})
            if j==0
                model=si.mi{i}.vesicleModel(:);
                j=1;
            else
                model(:,i)=si.mi{i}.vesicleModel(:);
            end;
        end;
    end;
    imacs(model);
    
    