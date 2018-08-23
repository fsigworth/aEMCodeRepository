% meFindDuplicateImages.m
% Search merged images for ones that show a large cross-correlation

% Have the user select some mi files: boilerplate
    [fname, pa]=uigetfile('*mi.*','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
    cd(rootPath);
nFiles=numel(fname);
%%
for i=1:nFiles
    mi=ReadMiFile([infoPath fname{i}]);
    m=BinImage(ReadEMFile(['Merged/' mi.baseFilename 'mz.tif']),4);
    if i==1
        ms=zeros([size(m) nFiles],'single');
        fms=zeros([size(m) nFiles],'single');
    end;
    ms(:,:,i)=m;
    fms(:,:,i)=fftn(m);
    imags(m);
    title(fname{i},'interpreter','none');
    drawnow;
end;
%%
maxsh=8;
ccs=zeros(nFiles,nFiles);
for i=1:nFiles
    for j=i+1:min(nFiles,i+maxsh)
        cc=GaussHP(fftshift(real(ifftn(fms(:,:,i).*conj(fms(:,:,j))))),.05);
        ccs(i,j)=max(cc(:));
        imags(cc);
        title(num2str([i j]));
        drawnow;
    end;
end;
imags(ccs);
%%
q=ccs(ccs>0);
figure(3);
hist(q,100);

sum(q>800) % CC larger than about 800 indicates a duplicated image.

save DuplicateData.mat ccs fname

%%
    infoPath1='Info/';
    for i=1:nFiles
        if any(ccs(i,:)>800)
            inds=find(ccs(i,:)>800);
%                 disp(num2str([i inds]));
                    origName=[infoPath1 fname{inds(1)}];
                    [pa,nm,ex]=fileparts(origName);
                    nm=[nm 'XD'];
                    newName=[AddSlash(pa) nm ex];
                    str=['!mv ' origName ' ' newName];
                    eval(str);
                    disp(str);
        end;
    end;           
            