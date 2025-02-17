function ShowSections( map, ctr, yAngle, options )
% function ShowSections( map, ctr, yAngle, options);
%  new version, previously called ShowSections2.
% Show three sections and projections of the 3d map.
% The optional ctr vector gives the [x y z] values at which sections are
% taken.  Default is FFT center, also taken when ctr=[].
% The optional yAngle is the angle in degrees of the y-projection (Y' axis)
% relative to the X-axis; default is 90.
% if map is nxnxnxnim, we show multiple lines, 6 images per line.
%
% Now handls complex maps too:
%   ShowSections(map,ctr,exponent)
% options.rowLabels doesn't work yet because mysubplot doesn't handle
% ylabel or xlabel correctly.

if nargin<4
    options=struct;
end;
    defOpts=struct;
    defOpts.showLines=1;
    defOpts.showLabels=0;
    defOpts.rowLabels={};
    defOpts.useMySubplot=1;
    defOpts.noAxis=0;
options=SetOptionValues(defOpts,options);




drawX=1;
drawY=1;

[nx, ny, nz, nim]=size(map);
if nargin<3
    yAngle=90;
else
    %     Rotations require size to be a multiple of 8.
    n=8*ceil(max([nx ny nz])/8);
    map=Crop(map,n);
    [nx, ny, nz, nim]=size(map);
end;
isStack=(nim>1);
isComplex=~all(imag(map(:))==0);
if isComplex
    if nargin<3
        exp=.5;
    else
        exp=yAngle(1);
    end;
    yAngle=90;
    SetComplex;
end;

ctr0=ceil(([nx ny nz]+1)/2);

if nargin<2 || numel(ctr)<3
    ctr=ctr0;
end;

showRowLabels=numel(options.rowLabels)>=nim;

xs=1:nx;
ys=1:ny;
zs=1:nz;

cx=ctr(1);
cy=ctr(2);
cz=ctr(3);

cx0=ctr0(1);
cy0=ctr0(2);
cz0=ctr0(3);

if yAngle>=90 || isComplex
    yxs=[cx cx];
    yxs0=[cx0 cx0];
    yys=[1 ny];
    yys0=yys;
    rotMap=map;
    yString='Y';
else
    tanY=tand(yAngle);
    yxs=[1 nx];
    yxs0=yxs;
    yys=cy+tanY*[1-cx nx-cx];
    yys0=cy0+tanY*[1-cx0 nx-cx0];
    m1=reshape(map,nx,ny,nz*nim);
    rotMap=rsRotateImage(m1,yAngle-90);
    rotMap=reshape(rotMap,nx,ny,nz,nim);
    yString='Y''';
end;

nc0=6;
if nim==1
    nr=3;
    nc=3;
else
    nr=nim;
    nc=6;
end;

for iSlice=1:nim
    % 1st column: X-Y plane
    subplot1(nr,nc,iSlice*nc0-5);
    imdraw(xs,ys,map(:,:,cz,iSlice));
    if options.showLines
        hold on;
        if drawX
            plot([1 nx ],[cy cy],'b-');
        end;
        % plot([cx cx],[1 ny],'g-');
        if drawY
            plot(yxs,yys,'g-');
        end;
        hold off;
    end
if options.noAxis
    axis off equal
end;
    if showRowLabels
        title(options.rowLabels{iSlice});
    elseif iSlice==1
            title('Sections XY');
    end;
    
    subplot1(nr,nc,iSlice*nc0-2);
    imdraw(xs,ys,sum(map(:,:,:,iSlice),3));
    if options.showLines
        hold on;
        %     plot([1 nx ],[cy0 cy0],'b-');
        %     % plot([cx0 cx0],[1 ny],'g-');
        %     plot(yxs0,yys0,'g-');
        if drawX
            plot([1 nx ],[cy cy],'b-');
        end;    % plot([cx cx],[1 ny],'g-');
        if drawY
            plot(yxs,yys,'g-');
        end;
        hold off;
    end;
    if options.showLabels
        xlabel('X');
        ylabel('Y');
    else
if options.noAxis
    axis off equal
end;
        if iSlice==1
            title('Projections XY');
        end;
    end;
    
    % X-Z plane
    subplot1(nr,nc,iSlice*nc0-4);
    imdraw(xs,zs,map(:,cy,:,iSlice));
    if options.showLines
        hold on;
        plot([1 nx ],[cz cz],'b-');
        plot([cx cx],[1 nz],'r-');
        hold off;
    end;
    if options.showLabels
        xlabel('X');
        ylabel('Z');
        title('Sections');
    else
if options.noAxis
    axis off equal
end;
        if iSlice==1
            title('XZ');
        end;
    end;
    
    subplot1(nr,nc,iSlice*nc0-1);
    imdraw(xs,ys,sum(map(:,:,:,iSlice),2));
    if options.showLines
        hold on;
        %     plot([1 nx ],[cy0 cy0],'b-');
        %     plot([cx0 cx0],[1 nz],'r-');
        if drawX
            plot([1 nx ],[cz cz],'b-');
        end;
        plot([cx cx],[1 nz],'r-');
        
        hold off;
    end;
    if options.showLabels
        xlabel('X');
        ylabel('Z');
        title('Projections');
    else
if options.noAxis
    axis off equal
end;
        if iSlice==1
            title('XZ');
        end;
    end;
    
    % Y-Z plane
    subplot1(nr,nc,iSlice*nc0-3);
    imdraw(ys,zs,rotMap(cx,:,:,iSlice));
    if options.showLines
        hold on;
        plot([1 ny ],[cz cz],'g-');
        plot([cy cy],[1 nz],'r-');
        hold off;
    end;
    if options.showLabels
        xlabel(yString);
        ylabel('Z');
    else
if options.noAxis
    axis off equal
end;
        if iSlice==1
            title('YZ');
        end;
        
    end;
    
    subplot1(nr,nc,iSlice*nc0);
    imdraw(xs,ys,sum(rotMap(:,:,:,iSlice),1));
    if options.showLines
        hold on;
        %     plot([1 ny ],[cz0 cz0],'g-');
        %     plot([cy0 cy0],[1 nz],'r-');
        plot([1 ny ],[cz cz],'g-');
        plot([cy cy],[1 nz],'r-');
        hold off;
    end;
    if options.showLabels
        xlabel(yString);
        ylabel('Z');
    else
if options.noAxis
    axis off equal
end;
        if iSlice==1
            title('YZ');
        end;
    end;
    if ~isStack
        %         rotMap=abs(rotMap);
        subplot(3,3,7);  % Plot 1D sections along x, y and z
        %         plot(xs,rotMap(:,cy,cz)','b-',...
        plot(xs,map(:,cy,cz)','b-',...
            ys,rotMap(cx,:,cz),'g-',....
            zs,squeeze(rotMap(cx,cy,:)),'r-','linewidth',2);
        axis([1 nx -inf inf]);
        legend('X','Y','Z');
        title('Line sections');
        %         set(gca,'color','k');
        
        subplot(3,3,8);  % Plot corresponding 1D projections
        plot(xs,squeeze(sum(sum(rotMap,3),2)),'b-',...
            ys,squeeze(sum(sum(rotMap,3),1)),'g-',....
            zs,squeeze(sum(sum(rotMap,2),1)),'r-','linewidth',2);
        axis([1 nx -inf inf]);
        title('Line projections');
        %         set(gca,'color','k');
    end;
end;

    function imdraw(xq,yq,m)
        if isComplex
            imacx(m,exp);
        else
            imags(xq,yq,m);
        end;
        
    end
    function subplot1(nr,nc,ind)
       mysubplot(nr,nc,ind,0,0,0,0);
%         subplot(nr,nc,ind);
    end;
    
end
