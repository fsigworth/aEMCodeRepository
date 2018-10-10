% rsCountVesicles
% Count all the picks in a directory full of info files
startEntry=1;
maxEntries=inf;
stride=1;
txtSize=12;
requireParticles=1;
d=dir('Info/');
imageFileSuffix='m.mrc';
vesicleFileSuffix='mv.mrc';
nPicks=0;
nEntries=min(numel(d),maxEntries);
ds=4;  % assume we're working with downsampled images
nmi=0;
miPicks=zeros(nEntries,1);
miDef=zeros(nEntries,1);
vesR=cell(nEntries,1);
vesS=cell(nEntries,1);
vesOk=cell(nEntries,1);

for i=startEntry:nEntries
    name=['Info/' d(i).name];
    
    if numel(name)<6 || ~strcmp(name(end-5:end),'mi.txt')
        continue;
    end;
    nmi=nmi+1;
    mi=ReadMiFile(name);
    miDef(nmi)=mi.ctf(1).defocus;
    msName=[mi.procPath mi.baseFilename imageFileSuffix];
    [msName,ok]=CheckForImageOrZTiff(msName);
    mvsName=[mi.procPath mi.baseFilename vesicleFileSuffix];
    [mvsName,ok2]=CheckForImageOrZTiff(mvsName);
    if ~(exist(msName,'file') && exist(mvsName,'file'     ))
        disp('no file');
        continue;
    end;
    m=ReadEMFile(msName);
    %     ds=mi.imageSize(1)/size(m,1);
    n=mi.imageSize(1)/ds;
    m=Downsample(m,n);
    mv=ReadEMFile(mvsName);
    mv=Downsample(mv,n);
    v=m-mv;
    mysubplot(131);
    imags(m);
    mysubplot(133);
    cla;
    mysubplot(132);
    imags(v);
    marksPresent=false;
    if isfield(mi.vesicle,'x') && numel(mi.vesicle.x)>0 % we have vesicles
        nv=numel(mi.vesicle.x);
        disp([name '  ' num2str(nv) ' vesicles.']);
        xs=double(mi.vesicle.x/ds+1);
        ys=mi.vesicle.y/ds+1;
        yls=double(ys+mi.vesicle.r(:,1)/ds);
        effAmps=zeros(nv,1);
        mi.vesicle.ok(:,5)=true;
        hold on;
        for j=1:nv
            text(xs(j),yls(j),num2str(j),'color','g',...
                'HorizontalAlignment','center');
        end;
        hold off;
        md=m;
        md(:,:,2)=v;
        k=1;
        b='R'; %%%%% s
        %             disp([(1:nv)' mi.vesicle.s(:,1,1)*1000 abs(mi.vesicle.r)]);
        %             disp(' ');
        % ---------------------interactive loop---------------------
        while b~='q' && b~='Q'
            %             disp([x y b]);
            switch b
                case 1
                    dists=hypot(x-xs,y-ys);
                    [minDist, k]=min(dists);
                    marksPresent=true;
                case {'n' 'R'}
                    mi.vesicle.ok(k,5)=true;
                    k=k+1;
                    marksPresent=true;
                case 'x'
                    mi.vesicle.ok(k,5)=false;
                    k=k+1;
                case 'N'
                    k=k-1;
            end;
            if k>nv
                b='q';
                k=nv;
            end;
            
            if k<1
                beep;
                k=1;
            end;
            %
            %             v1=meMakeModelVesicles(mi,960,k);
            effAmps(k)=1e4*EstimateImageAmplitude(mi,k,ds);
            v1=0;
            md(:,:,1)=m-v1;
            for iw=1:2
                mysubplot(1,3,iw);
                
                imags(md(:,:,iw));
                hold on;
                if iw==2
                    for j=1:nv
                        if mi.vesicle.ok(j,5)
                            text(xs(j),yls(j),num2str(j),'color','g',...
                                'HorizontalAlignment','center','FontSize',txtSize);
                        else
                            text(xs(j),yls(j),num2str(j),'color','b',...
                                'HorizontalAlignment','center','FontSize',txtSize);
                        end;
                        
                        text(xs(k),yls(k),num2str(k),'color','r', ...
                            'HorizontalAlignment','center','FontSize',txtSize);
                    end;
                end;
                r1=mi.vesicle.r(k,:)/ds;
                [xc,yc]=CircleLineSegments(mi.vesicle.r(k,:)/ds,min(10,100/r1(1)));
                xc=double(xc+mi.vesicle.x(k)/ds+1);
                yc=double(yc+mi.vesicle.y(k)/ds+1);
                plot(xc,yc,'r-','HitTest','off');
                hold off;
                
                
            end;
            
            subplot(1,3,3);
            hold off;
            hold on;
            oks=logical(mi.vesicle.ok(:,5));
            yVals=effAmps;
            rs=double(mi.vesicle.r(:,1));
            mxr=max(rs);
            rs1=rs+.01*mxr;
            %                 yVals=mi.vesicle.s(:,1,1);
            yVals(k+1:end)=NaN;
            for j=1:k
                if oks(j)
                    plot(rs(j),yVals(j),'bo');
                    text(rs1(j),yVals(j),num2str(j),'color','b','fontsize',txtSize);
                else
                    plot(rs(j),yVals(j),'ro');
                    text(rs1(j),yVals(j),num2str(j),'color','r','fontsize',txtSize);
                end;
            end;
            %                 plot(mi.vesicle.r(k,1),yVals(k),'k+');
            hold off;
            axis([0 max(rs) 0 max(effAmps)]);
            grid on;
            
            
            %                         hold off;
            txt=['(' num2str(k,3) ')  a=' num2str(effAmps(k),3) '  r=[ ' num2str(abs(mi.vesicle.r(k,[1 3:end])),3) ' ] '];
            title(txt);
            disp(txt);
            %             disp(b)
            if b~='R'
                [x,y,b]=Myginput(1);
            else
                drawnow;
            end;
        end; % while b
        % ------------------------------------------
        %%
        oks=effAmps>20 ;
        disp('Making model vesicles...');
        v1=meMakeModelVesicles(mi,960,find(oks));
        disp('done.');
        mysubplot(131);
        imags(m-v1);
        mysubplot(132);
        imags(v1);
        hold on;
        for j=1:nv
            if oks(j)
                text(xs(j),yls(j),num2str(j),'color','b', ...
                    'HorizontalAlignment','center','FontSize',txtSize);
                [xc,yc]=CircleLineSegments(mi.vesicle.r(j,:)/ds,min(10,100/r1(1)));
                xc=double(xc+mi.vesicle.x(j)/ds+1);
                yc=double(yc+mi.vesicle.y(j)/ds+1);
                plot(xc,yc,'b-','HitTest','off');
            end;
        end;
        hold off;
        
        
        
        subplot(133);
        cla;
        hold on;
        for j=1:nv
            if oks(j)
                plot(rs(j),yVals(j),'bo');
                text(rs1(j),yVals(j),num2str(j),'color','b','fontsize',txtSize);
            else
                plot(rs(j),yVals(j),'ro');
                text(rs1(j),yVals(j),num2str(j),'color','r','fontsize',txtSize);
            end;
        end;
        hold off;
        grid on;
        
        mysubplot(131)
        imags(m-v1);
        hold on;
        for j=1:nv
            if oks(j)
                %                 text(xs(j),yls(j),num2str(j),'color','b',...
                %                     'HorizontalAlignment','center','FontSize',txtSize);
                %                 text(xs(j),yls(j),num2str(j),'color','b', ...
                %                     'HorizontalAlignment','center','FontSize',txtSize);
            else
                %                 text(xs(j),yls(j),num2str(j),'color','r',...
                %                     'HorizontalAlignment','center','FontSize',txtSize);
                text(xs(j),yls(j),num2str(j),'color','r', ...
                    'HorizontalAlignment','center','FontSize',txtSize);
                [xc,yc]=CircleLineSegments(mi.vesicle.r(j,:)/ds,min(10,100/r1(1)));
                xc=double(xc+mi.vesicle.x(j)/ds+1);
                yc=double(yc+mi.vesicle.y(j)/ds+1);
                plot(xc,yc,'r-','HitTest','off');
                
            end;
        end;
        hold off;
        medianVesAmp=median(effAmps(oks))
        %         pause
        
        %%
        
        %         disp(' ');
        
        
        %         if size(mi.particle.picks,1)>0
        %             flags=mi.particle.picks(:,3);
        %             num=sum(flags>=16 & flags<48);
        %             nPicks=nPicks+num;
        %             miPicks(nmi)=num;
        %         end;
        %         vesR{nmi}=mi.vesicle.r;
        %         vesS{nmi}=mi.vesicle.s;
        %         vesOk{nmi}=mi.vesicle.ok & ~(requireParticles & size(mi.particle.picks,1)<1);
        %         if mod(nmi,stride)==0 || i>nEntries-3  % print out every 'stride' entries
        %             disp([num2str(nmi) ' ' d(i).name '  ' num2str([num nPicks])]);
        %             disp(abs(vesR{nmi}(1,:)));
        %         end;
        %
    end; % if numel(mi.vesicle.x)
    disp(' ');
    if b=='Q'
        break;
    end;
end; % for i
% miPicks=miPicks(1:nmi);
% miDef=miDef(1:nmi);
% subplot(313)
% hist(miPicks,100);
% xlabel('Particles per micrograph');
% ylabel('Frequency');
% subplot(312)
% plot(miPicks);
% ylabel('Particles per micrograph');
% xlabel('Micrograph');
% subplot(311);
% plot(miDef);
% ylabel('Defocus, \mum');
% xlabel('Micrograph');
disp('Done.');


return

% %%
% % Estimate vesicle amplitude from image integral
% k=47;
% ds=4;
% ves=meMakeModelVesicles(mi,mi.imageSize(1)/ds,k,0,0);
% imags(ves);
% sv=sum(-ves(:))
% a0=sum(mi.vesicleModel);
% a1=sqrt(2*pi)*mi.vesicle.extraSD*(mi.vesicle.extraPeaks*0+1);
% as=[a0 a1];
% s=squeeze(mi.vesicle.s(k,1,:));
% estA=4*pi*(mi.vesicle.r(k,1))^2*(as*s)/ds^2
%
% effVesAmp=ds^2*estA/(4*pi*(mi.vesicle.r(k,1)^2)*a0)
% effVesAmp0=ds^2*sv/(4*pi*(mi.vesicle.r(k,1)^2)*a0)


function effAmp=EstimateImageAmplitude(mi,k,ds)
% Estimate image amplitude from normalized image
ves=meMakeModelVesicles(mi,mi.imageSize(1)/ds,k,0,0);
sv=sum(ves(:));
mi1=mi;
mi1.vesicle.s(k,:,:)=0;
mi1.vesicle.s(k,1,1)=1;
ves1=meMakeModelVesicles(mi1,mi.imageSize(1)/ds,k,0,0);
sv1=sum(ves1(:));
effAmp=sv/sv1;
end

% %% Compute coherence function
% vpc0=meMakeModelVesicles(mi,960,0,1,1);
% v1=mv+vpc0;
% mi1=mi;
% mi1.vesicle.s(:,:,2:3)=0; % delete the extra peaks
% vp0=meMakeModelVesicles(mi1,960,28,0,0);
% fv=fftshift(fftn(v));
% fv0=fftshift(fftn(vp0));
% %%
% k=1e6;
% H=(fv.*fv0)./(abs(fv0).^2+k);
% figure(2);
% SetComplex;
% imacx(GaussFilt(H,.05),.3);
