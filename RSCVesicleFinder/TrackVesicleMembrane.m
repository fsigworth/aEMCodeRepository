function mi2=TrackVesicleMembrane(h);
% function mi2=TrackVesicleMembrane(h);
% Given the handle from Vesicle_finding_GUI or VesicleFinder, use the initial
% round vesicle models and try to track the membrane from the vesicle model in
% polar coordinates, filling all high-order terms in mi2.vesicle.r The returned
% structure mi2 has updated vesicle radius fields and the 1st element of the
% amplitude, i.e. mi.vesicle.s(ind,1).  Higher-order terms of s remain zero.
% The function reads the following fields from the structure h:
    % mi
    % filtImage
    % maskIndex
    % oldFilterFreqs
    % h.sav.vesicleRadii
    % h.sav.vesicleAmps
% The returned mi2 structure has updated values exclusively in mi2.vesicle

displayOn=0; % Show intermediate results on a second figure
maxSlopeSq=1;
maxClosure=5;
maxS1Ratio=10;
minRadiusA=80; % vesicles smaller than this won't be tracked at all.
ampRefitScale=0.75;  % Re-fitted amplitudes are high by about 30%
startInd=1;     % Which vesicle index to start with (for debugging)
minRDecay=.2;

% ds1=1;  % no further downsampling
mi=h.mi;
mi2=mi;  % output mi structure

m=h.filtImage;
n=size(m);
msk=meGetMask(mi,n,1:h.maskIndex);
m=m.*msk;
% m=Downsample(m,size(m)/ds1);
n=size(m);

effCT=meGetEffectiveCTF(mi,n,h.ds0);
filtFreqs=h.oldFilterFreqs;

if displayOn
    oldFig=gcf;
    figure(2);
end;
    ds=h.ds0;  % Net downsampling of our working image of size n
    imctr=ceil((mi.imageSize+1)/2); % center of original micrograph
    pixA=mi.pixA*ds;
% compute the effective CTF along with filter functions
    effCTF=effCT.*GaussHPKernel(n,filtFreqs(1)*pixA)...
        .*Gaussian(n,2,sqrt(1/(2*log(2)))*filtFreqs(2)*n);

    scl.n=n;
    scl.ds=ds;
    scl.dsShift=h.ds0Shift;

nv=numel(mi.vesicle.x);
numRefined=0;

for ind=startInd:nv;
%   Pick up the vesicle x,y and radius in our image
    cx=(mi.vesicle.x(ind)+h.ds0Shift(1))/ds+1;
    cy=(mi.vesicle.y(ind)+h.ds0Shift(2))/ds+1;
    cr=mi.vesicle.r(ind,1)/ds;
    if ~mi.vesicle.ok(ind,1) || cr<minRadiusA/pixA % skip nonexistent or small vesicles
        continue
    end;
    % r0=ceil(shiftWidthA/pixA); % shift of polar coordinates to include the peak
%   Make a centered vesicle model
    mi1=mi;
    mi1.vesicle.x=imctr(1);
    mi1.vesicle.y=imctr(2);
    mi1.vesicle.r=mi.vesicle.r(ind,1);
    vModel=meMakeModelVesicles(mi1,scl,1,0,0);
    vModel=real(ifftn(fftn(vModel).*ifftshift(effCTF))); % apply ctf
    
    % Extract the vesicle from the image, and the model
    r0=round(cr);
    nx=NextNiceNumber(8*r0);  % make a box with room for 4 x nominal radius
    mx=ExtractImage(m,round([cx cy]),nx);
    vx=ExtractImage(vModel,ceil((n+1)/2),nx);
    [x,y]=CircleLineSegments(mi1.vesicle.r/ds,10);
    xi=double(round(x+nx/2+1));
    yi=double(round(y+nx/2+1));
    if displayOn
        subplot(231);
        imags(mx);
        hold on;  % draw the initial vesicle path
        plot(xi,yi,'b-');
        hold off;
        title(ind);
    end;
    maskOverlap=any(xi<1|xi>nx|yi<1|yi>nx);
    if ~maskOverlap
        for i=1:numel(xi)
            maskOverlap = maskOverlap || mx(xi(i),yi(i))==0;
        end;
    end;
        
    mP=ToPolar(mx);
    vP=ToPolar(vx);
%     Make room for decreases in radius by shifting the radial coordinate
    ccP=circshift(real(ifft(fft(mP).*conj(fft(vP)))),r0);
    nt=2*nx;
    nr=nx/2;
    
    % Find the origin for theta: pick the largest region near r0 value.
    [mxv,mxi]=max(ccP(1:2*r0,:)); % get the maximum at each theta
    mxCount=(mxi>0.8*r0 & mxi<1.2*r0)';
    [mxv,t0]=max(GaussFilt(mxCount,1/(2*nx)));  % t0 is the theta origin of our search.
    ccP0=circshift(ccP,[0 -t0]);
    
    % estimate noise
    sds=std(ccP0(round(r0+0.5*cr):round(r0+1.5*cr),:));
    sigma=median(sds);  % very rough estimate of noise
    
    % field of observable probabilites
    pccP0=exp(ccP0/sigma);
    % for i=1:nx*2
    %     pccP0(:,i)=pccP0(:,i)/sum(pccP0(:,i));
    % end;
    
    % initial probability at the theta origin
%     p0=circshift(Gaussian(nx/2,1,r0/2,nx/4+1),round(r0-nx/4-1));
    p0=circshift(Gaussian(nx/2,1,r0/4,nx/4+1),round(r0-nx/4));
    
    %% Track the path
    a=.1;
    A=[a 1-2*a a]; % transition matrix
    path=ViterbiPath(p0,pccP0,A);
    
    path=path-1;  % shrink by 1 pixel
    
    if displayOn
        subplot(234);
        imags(ccP0');
        hold on;
        plot(path,'linewidth',2);
        hold off;
        title('pccP0');
    end;
    q=diff(diff(GaussFilt(path,.02)));
    slopesq=max(abs(q));
    closure=path(end)-path(1);
    % disp([slopesq*10 closure]);
    if displayOn
        subplot(236);
        plot(q);
    end;
    pathsh=circshift(mod(path+cr-r0-1,nx/2)+1,t0);

%     Make a Fourier expansion of the path, shifting the origin to force
%     the 1st-order term to zero.
    nTerms=4;
    thetas=(0:nt-1)'*2*pi/nt;
    px=cos(thetas).*pathsh;
    py=sin(thetas).*pathsh;
    % plot(px,py);
    % Get the new center point
    cx1=(max(px)+min(px))/2;
    cy1=(max(py)+min(py))/2;
    
    
    for j=1:5 % iterations of origin shifting
        th1=atan2(py-cy1,px-cx1);  % thetas of our points
        th1(th1<0)=2*pi+th1(th1<0);
        th1(1)=0;
        r1=hypot(py-cy1,px-cx1);
        % Get the Fourier expansion about this center
        f=zeros(nt,2*nTerms+1);
        f(:,1)=1;
        for i=2:2:2*nTerms
            w=i/2; % [0 1 1 2 2 3 3...]
            f(:,i)=cos(w*th1);
            f(:,i+1)=sin(w*th1);
        end;
        warning('off');
        a=LinLeastSquares(f,r1);
        warning('on');
        cx1=cx1+a(2); % cos(theta) coeff gives x shift
        cy1=cy1+a(3); % sin(theta) coeff gives y shift
        %     disp([j a(2) a(3)]);
        if hypot(a(2),a(3))<.1
            break
        end;
        %     a(2:3)=0;
        %     plot(th1, [f*a r1]);
        %      plot([px-cx1 f*a.*cos(th1)],[py-cy1 f*a.*sin(th1)]);
        %     pause(0.2)
    end;
    if displayOn
        subplot(235);
        plot([px-cx1 f*a.*cos(th1)],[py-cy1 f*a.*sin(th1)]);
    end;
    % r=[a(1) 1*a(4:2:end)'-1*1i*a(5:2:end)'];
    r=[a(1) 1*a(2:2:end)'-1*1i*a(3:2:end)'];
    
    slopeSqVal=100*slopesq/sqrt(a(1));
    if displayOn
        subplot(236);
        title(['Slope ' num2str(100*slopesq/sqrt(a(1)))]);
%         title(num2str([a(1) 100*slopesq/sqrt(a(1)) closure]));
    end;
    
    % disp(r);
    mi1.vesicle.r=r*ds;
    mi1.vesicle.s=mi.vesicle.s(ind,1);
%     mi1.vesicle.s=.001;
    mi1.vesicle.s(1,nTerms)=0; % pad with zeros
    % mi1.vesicle.extraS=mi1.vesicle.s;
    mi1.vesicle.ok=true(1,4);
    v1=meMakeModelVesicles(mi1,scl,1,0,0);
    v1ctr=round([(mi1.vesicle.x(1)+h.ds0Shift(1))/ds+1 (mi1.vesicle.y(1)+h.ds0Shift(2))/ds+1]);
    v1Model=real(ifftn(fftn(v1).*ifftshift(effCTF))); % apply ctf
    
    if displayOn
        v1x=ExtractImage(v1Model,v1ctr,nx);
        
        subplot(232);
        imags(v1x);
        hold on;
        plot(px-cx1+nx/2,py-cy1+nx/2,'y-');
        hold off;
    end;
    
%     Modify the amplitude to fit the vesicle image
    v2x=ExtractImage(v1Model,round(v1ctr-[cx1 cy1]),nx);
    s1=ampRefitScale*v2x(:)'*mx(:)/(v2x(:)'*v2x(:));
%     Explicit lin least squares yields essentially the same result:    
%     n2x=numel(v2x);
%     f=zeros(n2x,2);
%     f(:,1)=sum(v2x(:));
%     f(:,2)=v2x(:);
%     sVals=LinLeastSquares(f,mx(:));
% disp([ind s1 sVals(2)]);
%     s1=sVals(2);
    
    if displayOn
        title(['amps ' num2str(1000*mi.vesicle.s(ind,1)*[1 s1])]);
        subplot(233);
        imags(mx-s1*v2x)
    end;
    
    % Check to see if criteria are met
    if ~maskOverlap && slopeSqVal<maxSlopeSq && closure<maxClosure && ~isnan(s1) && s1>1/maxS1Ratio && s1<maxS1Ratio
        numRefined=numRefined+1;
        ner=numel(mi1.vesicle.r);
        mi2.vesicle.r(ind,1:ner)=mi1.vesicle.r;
        nes=numel(mi1.vesicle.s);
        mi2.vesicle.s(ind,1:nes)=s1*mi1.vesicle.s;
        mi2.vesicle.x(ind)=mi.vesicle.x(ind)+ds*cx1;
        mi2.vesicle.y(ind)=mi.vesicle.y(ind)+ds*cy1;
        % Mark our vesicle ok if it meets the criteria.  ok(ind,3) is set
        % if tracking was successful.
        r01=mi2.vesicle.r(ind,1)*mi2.pixA;
        s01=mi2.vesicle.s(ind,1);
rDecay=minRDecay+(1-minRDecay)/(1+(r01/200)^2); % help for big vesicles.

        flag=(r01>=h.sav.vesicleRadii(1) && r01<=h.sav.vesicleRadii(2) ...
            && s01>=h.sav.vesicleAmps(1)*rDecay && s01<=h.sav.vesicleAmps(2));
        mi2.vesicle.ok(ind,:)=[1 flag 1 0];

        % blank all vesicles whose centers overlap with our new vesicle.
        vesMask=VesicleMaskGeneral(n,mi2.vesicle.r(ind,:)/ds,0,[mi2.vesicle.x(ind) mi2.vesicle.y(ind)]/ds+1);
        nBlanked=0;
        for i=1:nv
            if i~=ind
                xi=max(1,min(n(1),round(mi2.vesicle.x(i)/ds)));
                yi=max(1,min(n(2),round(mi2.vesicle.y(i)/ds)));
                if vesMask(xi,yi)
                    mi2.vesicle.ok(i,:)=false;
                    nBlanked=nBlanked+1;
                end;
            end;
        end;
        title(['updated with ' num2str(nBlanked) ' blanked.']);
    else
        if displayOn
            subplot(235);
            title('----');
            subplot(234);
            hold on;
            plot(path,'r-','linewidth',2);
            hold off;
        end;
    end;
    if displayOn
        drawnow;
    end;
end; % loop over vesicle index

if displayOn
    figure(oldFig);
end;


disp([num2str(numRefined) ' vesicles refined of ' num2str(nv)]);

return





% Compute center of mass
ts=1:nt;
expT=ts*path/sum(path);
expR=sum(path)/nt;

for i=1:2*nx
    bImage(round(pathsh(i)),i)=1;
end;
%
imags(ccP0');
hold on;
plot(path)
hold off;

bImage=zeros(nx/2,2*nx,'single');
pathsh=circshift(path,t0);
for i=1:2*nx
    bImage(round(pathsh(i)+cr-r0),i)=1;
end;
% bImage(round(expR),1)=2;
imags(bImage)
pause(0.1)
bx=ToRect(bImage);
imags(bx+mx);

% Get the approximate center
bBits=circshift(bx>.1,[0 0]);
[X,Y]=ndgrid(-nx/2:nx/2-1);
dx=X(:)'*bx(:)/sum(bx(:))
dy=Y(:)'*bx(:)/sum(bx(:))
hold on;
ctx=ceil((nx+1)/2);
plot(dx+ctx,dy+ctx,'y+');
hold off;
% Convert to center again
vmP=ToPolar(circshift(bx,round([-dx -dy])));
imags(vmP);
[vals,mpath]=max(vmP);  % mPath is radius as a function of angle.
fpath=fft(mpath);
plot(abs(fpath(1:7)));
fp0=fpath;
fp0(2:7)=2*fp0(2:7);
fp0(8:end)=0;
smpath=real(ifft(fp0));
plot([mpath' smpath']);



