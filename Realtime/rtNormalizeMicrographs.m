% rtNormalizeMicrographs.m
% sitting in the Micrograph folder, make normalized images and write them into
% the Merged folder.

figure(2);
pwd
d=dir;
for i=1:numel(d)
    if ~d(i).isdir && strndcmp(d(i).name,'.mrc')
        nmi=d(i).name;
        nmo=nmi;
        % replace 'ali' with 'm'
        p=strfind(nmo,'ali');
        if numel(p)>0
            nmo(p(end):p(end)+2)=[];
        end;
        [pa,nm,ex]=fileparts(nmo);
        nmo=['../Merged/' nm 'm' ex];
        disp(nmo);
        [m,s]=ReadMRC(nmi);
%         s.pixA
        m=m/mean(m(:))-1;
        WriteMRC(m,s.pixA,nmo);
        imags(BinImage(m,4));
        title(nmo,'interpreter','none');
        drawnow;
%     return
    end;
end;