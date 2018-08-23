% MiDefocusTrend.m
% load AllMis.mat
MiLoadAll;
%%
clf
nmi=numel(mis);
stride=Step125(nmi/25);

ds=zeros(nmi,2);
for i=1:nmi
    mi=mis{i};
    if isfield(mi,'ctf') && isfield(mi.ctf,'defocus')
    for j=1:numel(mi.ctf)
        ds(i,j)=mi.ctf(j).defocus;
    end;
    end;
end;
plot(ds,'.-');
numVals=sum(ds(:,1)>0);
disp([num2str(numVals) ' nonzero defocus values.']);
for i=1:stride:nmi
    text(i,8,mis{i}.movieFilename{1}(1:end-4),'rotation',90,'fontsize',10,'interpreter','none');
end;

xlabel('Micrograph');
ylabel('Defocus, \mum');
grid on;
disp('defocus display done.');

return
%%
%fo=fopen('DefocusList.txt','w');
fo=1;
for i=1:nmi
    mi=mis{i};
    if isfield(mi,'ctf') && isfield(mi.ctf,'defocus')
        def=mi.ctf(1).defocus;
    else
        def=0;
    end;
    if ~isfield(mi,'movieFilename')
        mi.movieFilename={'**'};
    end;
    %fprintf(fo,'%4g  %8.4f  %s\n',i,def,mi.movieFilename{1});
end;
%fclose(fo);
