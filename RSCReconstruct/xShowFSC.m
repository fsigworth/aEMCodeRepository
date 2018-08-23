% xShowFSC.m
% Given a reconstruction directory, show the volume and fsc
% fscs=[];
%
mrcDir='mrc/';
jpegDir='jpeg/';

doPause=0;
showMasked=0;
doAlign=0;
startFig=5;
d=dir(mrcDir);
i=0;
as=zeros(0,4);
for ind=3:numel(d)
    [pa,nm,ex]=fileparts(d(ind).name);
    if strcmp(ex,'.mrc')
        a=sscanf(d(ind).name,'%c%d%c%c%d%s');  % e.g. i01av01.mrc
        if numel(a)==9
            i=i+1;
            as(i,:)=[a(2) a(3)-96 a(5) ind];
        end;
    end;
end;
if i==0
    error('No .mrc files found.');
end;
firstIter=min(as(:,1));
lastIter=max(as(:,1));
nTwins=max(as(:,2));
nVols=max(as(:,3));

% firstIter=lastIter;
lastIter
firstIter=MyInput('Starting iteration',firstIter);
firstIter=min(firstIter,lastIter);

for ind=firstIter:lastIter
    is=find(as(:,1)==ind);
    if numel(i)<1 break; end;
    %%
    for i=is'
        figure(startFig+i-is(1)+1);
        set(gcf,'menubar','none');
        pos=get(gcf,'position');
        pos(3:4)=[560 480];
        set(gcf,'position',pos);
        vi=i-is(1)+1;
        name=[mrcDir d(as(i,4)).name];
        [pa,nm,ex]=fileparts(name);
        disp(name);
        if vi==1
            [v,s]=ReadMRC(name);
            n=size(v,1);
            ctr=floor(n/2+1);
            volMask=min(1,fuzzymask(n,3,[.32 .32 .15]*n,.1*n,[ctr ctr ctr+.18*n])...
                +fuzzymask(n,3,0.22*n,0.1*n,[ctr ctr ctr-.06*n]));
        else
            v(:,:,:,vi)=ReadMRC(name);
        end;
        ShowSections(v(:,:,:,vi),[],45);
        title(name);
        drawnow;
        
        if showMasked && ind==lastIter  % write out the masked image
            WriteMRC(v(:,:,:,vi).*volMask,s.pixA,[nm 'msk' ex]);
            figure(i-is(1)+1);
            ShowSections(v(:,:,:,vi).*volMask,[],45);
            title(name);
        end;
        
        if as(i,2)==2 && i>1  % we have a twin
            vmsk=volMask;
            if doAlign
                [vol2ali,tz,gamma,mirror]=reAlignVolumes(v(:,:,:,vi-1),v(:,:,:,vi));
                disp([tz gamma mirror]);
                ShowSections(vol2ali,[],45);
                nm=[nm 'ali'];  % change the stored jpeg name                
                title(nm);
                WriteMRC(vol2ali,s.pixA,[mrcDir nm '.mrc']);
            else
                vol2ali=v(:,:,:,vi);
            end;
            n=size(v,1);
            fs=(0:n/2-1)/(n*s.pixA);
            fsc=FSCorr2(vmsk.*v(:,:,:,vi-1),vol2ali);
            mysubplot(339);
            plot(fs,[fsc fsc*0]);
            figure(startFig+i-is(1)+3);
            set(gcf,'menubar','none');
            pos=get(gcf,'position');
            pos(3:4)=[560 480];
            set(gcf,'position',pos);
            
            mergeV=v(:,:,:,vi-1)+vol2ali;
%            WriteMRC(mergeV.*volMask,s.pixA,[mrcDir nm '.mrc']);
            ShowSections(mergeV.*volMask,[],45);
        end;
        jpegName=[jpegDir nm '.jpg'];
        disp(['Writing ' jpegName]);
        print(jpegName,'-djpeg','-r150');
        
    end;
    
    
    
    if doPause
        pause;
    end;
end;


return


iTxt=['i' sprintf('%02d',ind)];
[v1,s]=ReadMRC(['mrc/' iTxt 'av01.mrc']);
%=======
for iter=50;
    for iVol=1:5
        iTxt=['i' sprintf('%02dav%02d.mrc',iter,iVol)];
        [v1,s]=ReadMRC(iTxt);
        figure(iVol);
        ShowSections(v1);
    end;
end;
%%




% >>>>>>> Stashed changes
%  v2=ReadMRC(['mrc/' iTxt 'bv01.mrc']);
%  [vol2ali,tz,gamma,mirror]=reAlignVolumes(v1,v2);
%  vol2ali=v2;
%
% disp(['iter ' sprintf('%3d',ind) '  alignment: ' sprintf('[%5.1f %5.1f  %1d ]',[tz gamma mirror])]);
%
%  n=size(v2,1);
%  cn=ceil(n/2+1);
% %   vmsk=fuzzymask(n,3,0.2*n,0.1*n,[cn cn cn-12]);
% %   vmsk=fuzzymask(n,3,0.3*n,0.1*n,[cn cn cn]);
% vmsk=min(1,fuzzymask(n,3,[.32 .32 .15]*n,.1*n,[cn cn cn+.18*n])+fuzzymask(n,3,0.22*n,0.1*n,[cn cn cn-.06*n])); figure(1);
%  ShowSections(vmsk.*(v1),[],45);
%  [p2,p1]=ParsePath(pwd);
%
%  figure(2);
%  ShowSections(vmsk.*vol2ali,[],45);
%
%  pause
% return
%
%
% load ri
%  subplot(332);
%  title(ri.siName,'interpreter','none');
%
% index=3;
% fscs(:,index)=fsc;
