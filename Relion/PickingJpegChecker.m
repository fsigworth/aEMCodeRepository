% PickingJpegChecker

jpegDir='Picker_jpegs/';
%  assume we've loaded fracs, to look at fractional overlaps.
% load allMis.mat

% % Intensities
% load Picking_8/allMis8_intens.mat
% nmi=numel(allMis);
% intens=zeros(10,1);
% for i=1:nmi
%     mi=allMis{i};
%     intens(:,i)=mi.ok(11:20);
% end;
% 

nmi=numel(allMis);
vals=zeros(nmi,15);
for i=1:nmi
    z=allMis{i}.ok(1:15);
    vals(i,:)=z;
end;
%  Set the selection criteria here:
q=1./(1+vals(:,1).^2); % switch function favoring overlap <=1
% vals(:,8) is the composite intensity
selVal=vals(:,8)+(8.4-vals(:,8)).*q/2;
selVal2=vals(:,1);
figure(2);
clf;
hist(selVal,200);

figure(3);
selVal2(selVal2>20)=NaN;
plot(selVal,selVal2,'.');

inds=find(selVal>8.2 & selVal<8.5 & defs<3.6 & ~isnan(selVal) & res<7);
% inds=find(selVal2<3);

figure(1);
clf;
set(1,'color',[.4 .4 .4]);
set(1,'menu','none');
set(1,'toolbar','none');
% inds=find(rej);
% inds=1:nmi;
nin=numel(inds);
disp([num2str(nin) ' indices']);
iptr=1;

while iptr<nin
    ind=inds(iptr);
    mi=allMis{ind};
    allMis{ind}.active=true;
    jpegName=[jpegDir mi.baseFilename '_i0.jpg'];
        tstr=sprintf('%04d %04d %6.3f %6.3f',iptr,ind,selVal(ind),selVal2(ind));
%         tstr=sprintf('%04d %6.3f',iptr,fracs(iptr));
%         tstr=sprintf('%04d file: %05d   def: %5.2f   res: %6.2f   FOM: %6.3f',...
%             iptr,ind,mi.ctf.defocus, mi.ctf.resLimit,mi.ctf.ccc);
    if exist(jpegName,'file')
        im=imread(jpegName);
        image(im);
%         fprintf('%d %g',ind,selVal(ind));
disp(tstr);n
    else
        iptr=iptr+1;
        continue;
        disp([num2str(ind) ' not found: ' jpegName]);
        imags(zeros(768));
    end
        title(tstr,'fontsize',16,'fontweight','normal','color','w');
    [x,y,b]=ginput(1);
    switch char(b)
        case 'i'
            tptr=MyInput('Pointer value ',iptr);
            iptr=tptr-1;
        case 'j' % junk
            allMis{ind}.active=false;
            fprintf(' X\n');
        case 'p' % go to previous
            iptr=iptr-2;
        case 'g' % go to a particular indey
            ind=MyInput('File index to load ',iptr);
            if ind>=1 && ind<nin
                iptr=ind-1;
            end;
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
 imags(zeros(768));    
 disp(' ');
 disp('Done.');
 return
 %%
 
 % Set the active flags
 for i=1:nmi
     mi=allMis{i};
     mi.active=~rej(i);
     allMis{i}=mi;
 end;
 return
 
%%  Set ok flags based on intensity
lowerLimit=8.1; % selects about half of the total. To get 1/3, lower limit=8.2.
upperLimit=8.6; 
for i=1:nmi
    allMis{i}.ok(2)=intens(8,i)>lowerLimit & intens(8,i)<upperLimit;
end;
save Picking_8/allMis8_intens2.mat allMis


 %%
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
