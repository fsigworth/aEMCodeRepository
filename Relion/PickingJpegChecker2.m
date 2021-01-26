% PickingJpegChecker2
% Combined with CountOverlappingVesicles

jpegDir='Picker_jpegs/';
%  assume we've loaded fracs, to look at fractional overlaps.
% load allMis.mat
nmi=numel(allMis);

figure(1);
set(1,'color',[.4 .4 .4]);
set(1,'menu','none');
set(1,'toolbar','none');
% inds=find(rej);
% inds=1:numel(fracs);
inds=1:numel(allMis);
nin=numel(inds);
disp([num2str(nin) ' indices']);
iptr=210;
b='0';
while iptr<nin
    
    ind=inds(iptr);
    mi=allMis{ind};
    allMis{ind}.active=true;
    jpegName=[jpegDir mi.baseFilename '_i0.jpg'];
    %         tstr=sprintf('%04d file: %05d   def: %5.2f   res: %6.2f   FOM: %6.3f',...
    %             iptr,ind,mi.ctf.defocus, mi.ctf.resLimit,mi.ctf.ccc);
    if exist(jpegName,'file')
        im=imread(jpegName);
        subplot(121);
        image(im);
    else
%         disp([num2str(ind) ' not found: ' jpegName]);
%         imags(zeros(768));
        tstr='';
    end;
    if ~isfield(mi,'vesicle') || numel(mi.vesicle.x)<1
            iptr=iptr+1;
        continue;
    end;
    subplot(122);
    fprintf('%d',ind);
    [ovImage,frac]=FractionOverlaps(allMis{iptr});
    imaga(128*ovImage);
    tstr=sprintf('%04d %4d',iptr,round(frac*100));
    allMis{ind}.ok(1)=100*frac;
    subplot(121);
    title(tstr,'fontsize',16,'fontweight','normal','color','w');
    drawnow;
    if b~='R'
        [x,y,b]=ginput(1);
    end;
    switch char(b)
        case 'i'
            tptr=MyInput('Pointer value ',iptr);
            iptr=tptr-1;
        case 'j' % junk
            allMis{ind}.active=false;
            fprintf(' X\n');
        case 'p' % go to previous
            iptr=iptr-2;
        case 'q' % quit (and set active=1...)
            break;
    end;
    if allMis{ind}.active
        fprintf('\n');
    end;
    iptr=iptr+1;
    if iptr<1
        beep
        ind=1;
    end;
end;

disp('Done.');
return
%%

% Set the active flags
for i=1:nmi
    mi=allMis{i};
    mi.active=~rej(i);
    allMis{i}=mi;
end;

% save the allMisSel.mat
disp('Saving allMis8_sel.mat');
save allMis8_sel.m allMis

%%
% save the individual mi files
infoDir='Info8_sel/';
CheckAndMakeDir(infoDir,1);
for i=1:nmi
    mi=allMis{i};
    miName=[infoDir mi.baseFilename 'mi.txt'];
    WriteMiFile(mi,miName);
    if mod(i,1000)==0
        disp([num2str(i) ' ' miName]);
    end;
end;

% store all the mi files
