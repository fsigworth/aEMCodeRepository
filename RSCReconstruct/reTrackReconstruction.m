% reTrackReconstruction.m
% Assuming that we are in a reconstruction directory, e.g.
% Reconstructions/Recon122n1/
% search for .mrc files, display them and compute the FSC, and store jpegs
% of the ShowSections figures.
restart=0;
startIter=1;

mrcDir='mrc/';
jpegDir='jpeg/';

doPause=0;
showMasked=0;
doAlign=0;

startFig=5;
figure(startFig);
pos0=get(gcf,'position');
for j=1:3
    figure(startFig+j-1);
    set(gcf,'menubar','none');
    pos=get(gcf,'position');
    pos(3:4)=[460 480];
    pos(1)=pos0(1)+500*(j-1);
    set(gcf,'position',pos);
end;
if ~restart
    as=zeros(0,4);
    oldAs=as;
    oldN=0;
    newInd=0;
    oldInd=0;
    flags=false(1,2);
end;

while 1
    %     disp('Reading dir');
    d=dir(mrcDir);
    for i=3:numel(d)
        [pa,nm,ex]=fileparts(d(i).name);
        if strcmp(ex,'.mrc')
            a=sscanf(d(i).name,'%c%d%c%c%d%s');  % e.g. i01av01.mrc
            if numel(a)==9
                a1=[a(2) a(3)-96 a(5)]; % [iter iTwin iVol]
                if ~any(all(repmat(a1,size(oldAs,1),1)==oldAs(:,1:3),2),1)  % unique new entry?
                    newInd=newInd+1;
                    as(newInd,:)=[a1 i];
                end;
            end;
        end;
    end;
    if newInd<=oldInd
        pause(10);
    else  % something new is found
        for ind=oldInd+1:newInd
            name=[mrcDir d(as(ind,4)).name];
            iter=as(ind,1);
            if iter<startIter
                disp(['skip ' name]);
            else
                iTwin=as(ind,2);
                iVol=as(ind,3);
                [pa,nm,ex]=fileparts(name);
                %          disp(name);
                iv=iTwin+iVol-1;  % figure index
                [v1,s]=ReadMRC(name);
                n=size(v1,1);
                if n~=oldN
                    flags=false(1,2);
                    oldN=n;
                    v=zeros(n,n,n,iv);
                    ctr=floor(n/2+1);
                    volMask=min(1,fuzzymask(n,3,[.32 .32 .15]*n,.1*n,[ctr ctr ctr+.18*n])...
                        +fuzzymask(n,3,0.22*n,0.1*n,[ctr ctr ctr-.06*n]));
                end;
                %             if iTwin>2
                %                 break  % skip 'c' volumes.
                %             end;
                flags(iTwin)=true;
                v(:,:,:,iv)=v1;
                if iTwin<3
                    figure(startFig+iv-1);
                    ShowSections(v(:,:,:,iv),[],45);
                    title(name);
                    drawnow;
                    
                    % Compute FSCs
                    if iTwin==2 && all(flags(1:2))  % we have a twin; compute and display FSC
                        vmsk=volMask;
                        n=size(v,1);
                        fs=(0:n/2-1)/(n*s.pixA);
                        fsc=FSCorr2(vmsk.*v(:,:,:,1),vmsk.*v(:,:,:,2));
                        [v2ali,tz,gamma,mirror]=reAlignVolumes(v(:,:,:,1),v(:,:,:,2).*vmsk);
                        disp([tz gamma mirror]);
                        fscAli=FSCorr2(v(:,:,:,1).*vmsk,v2ali);
                        subplot(339);
                        plot(fs,[fsc fscAli fsc*0]);
                        title(num2str([tz gamma mirror]));
                        aliPresent=1;
                    else
                        aliPresent=0;
                    end;
                    jpegName=[jpegDir nm '.jpg'];
                    disp(['Writing ' jpegName]);
                    print(jpegName,'-djpeg','-r150');
                    if aliPresent  % show another figure
                        figure(startFig+iv);
                        ShowSections(v2ali+v(:,:,:,iv-1),[],45);
                        subplot(339);
                        plot(fs,[fsc fscAli fsc*0]);
                        title(num2str([tz gamma mirror]));
                        jpegName=[jpegDir nm 'AliC.jpg'];
                        disp(['Writing ' jpegName]);
                        print(jpegName,'-djpeg','-r150');
                    end;
                end;
            end;
        end;
    end;
    oldAs=as;
    oldInd=newInd;
end % while


