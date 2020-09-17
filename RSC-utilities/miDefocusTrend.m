% MiDefocusTrend.m
load AllMis.mat
% MiLoadAll;
%%
mis=allMis;

clf
nmi=numel(mis);
stride=Step125(nmi/25);
labelY=0.1;

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
        imgName=mis{i}.baseFilename;
    text(i,labelY,imgName,'rotation',90,'fontsize',10,'interpreter','none');
end;

xlabel('Micrograph');
ylabel('Defocus, \mum');
grid on;
disp('defocus display done.');

return

%% Find the micrographs closest to a given defocus value.
targetDef=2.26;
numToChoose=10;
outDir1='Merged_sel/';
outDir0='Merged_sm_sel/';
outDir2='Info_sel/';
CheckAndMakeDir(outDir1);
defs=ds(:,1);
devs=abs(defs-targetDef);
[sdevs,sinds]=sort(devs);
cinds=sinds(1:numToChoose);
for i=1:numToChoose
    baseName=mis{i}.baseFilename;
    disp([num2str([i cinds(i)] ) '  '  num2str(defs(cinds(i) ),3) '  ' baseName]);
    m=ReadMRC(['Merged_sm/' baseName 'ms.mrc']);
    imags(GaussFilt(m,.1));
    pause (0.2);
    str0=['cp Merged_sm/' baseName 'ms.mrc ' outDir0];
    str1=['cp Merged/' baseName 'm.mrc ' outDir1];
    str2=['cp Info/' baseName 'mi* ' outDir2];
    disp(str0);
    system(str0);
    disp(str1);
    system(str1);
    disp(str2);
    system(str2);
end;




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
