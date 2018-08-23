function [p,shift,T,imtc]=StereoFitIteration(p,shift,im0,imt,rs)
%     rs is a struct of parameters
%         ds: overall downsampling
%         fc, f0, ft: filter freqs for cc, im0 and imt
    fm0=fftn(im0);
    nw=size(imt);
    ctr=ceil((nw+1)/2);
    
    theta=p(1);  % Tilt axis
    ang1=p(2);   % tilt angle
    T0=[1/cosd(ang1) 0 0; 0 1 0; 0 0 1];
    R0=[cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1];
    T=p(3)*R0*T0*inv(R0);  % apply this transform to tilted image
%     shift=max(min(shift,rs.mxShift),-rs.mxShift);
    T(1:2,3)=-shift;
    imtc=AffineTransform(imt,T);  % expand it.
    %                 imtc=circshift(imtc,-round(shift.*nw));
    im0f=GaussFilt(im0,rs.f0*rs.ds);
    subplot(224);
    imacs(rot90(im0f));
    drawnow;
    
    cc=(fftshift(real(ifftn(fftn(imtc).*conj(fm0))))/prod(nw));

    if rs.lpf<1
        cc=GaussFilt(cc,rs.lpf);
    end;
    cc=rs.mask.*(GaussHP(cc,rs.hpf*rs.ds)-rs.prior);
    
    %         Show the cross-correlation, and mark the maximum
    crcc=cc;
    subplot(223);
    scl=rs.pixA0*rs.ds/10;
    xs=(1:nw)*scl;
    imacs(xs,xs,rot90(crcc));
    hold on;
    [pk, ixc, iyc]=max2di(cc);
    ix=ixc-ctr(1);
    iy=iyc-ctr(2);
    plot(-iyc*scl,ixc*scl,'b+');
    hold off;
    xlabel('Displacement, nm');
    % Find the overall maximum, and shift the tilted image correspondingly.
%     [pk, ix, iy]=max2di(cc);
    shift=shift+([ix iy]./nw);

    %             imts=circshift(imtc,-round([ix-ctr iy-ctr]));
    imtf=GaussFilt(imtc,rs.ft*rs.ds);
    subplot(222);
    imacs(rot90(imtf));
    title(sprintf('%8.3f  ',[p shift]));
    subplot(224);
    imacs(rot90(imtf));
    subplot(221);
    imacs(rot90(im0f));
    drawnow;
      p=Simplex(-pk);

