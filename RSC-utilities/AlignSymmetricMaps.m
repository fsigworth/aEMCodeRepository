function [rot, zsh, m1dr]=AlignSymmetricMaps(m1,m2,ds,rot,zsh);
% function [rot, zsh]=AlignSymmetricMaps(m1,m2,ds,rot,zsh);
% rotate m1 about the z-axis, and shift along the z-axis, to maximize cross
% correlation, and return the angle in degrees and shift in pixels. ds is the downsampling factor to speed things up. If ds is a
% vector, we repeat the alignment with each given factor; default is
% [4 2 1]. rot and zSh are initial values for rotation (ccw in degrees) and
% z-shift in original pixels; defaults are 0 and 0. The returned map m1dr
% is the optimized transformation of m1.

if nargin<3 || numel(ds)==0
    ds=[4 2 1];
end;
if nargin<4
    rot=0;
end;
if nargin<5
    zsh=0;
end;

n0=size(m1,1);
numIters=100;
minAngle=.001;
displayIters=1;

for j=1:numel(ds)
    ds1=ds(j);
    n1=round(n0/ds1); % hoping n0/ds1 is integer!
%     ctr=ceil((n1+1)/2);
    m1d=Downsample(m1,n1);
    m2d=Downsample(m2,n1);

    rotStep=ds1;
    zsh=zsh/ds1;

    P=Simplex('init',[zsh rot],rotStep);
        r=P(2);
        zsh=P(1);
    for i=1:numIters
%         cc=GaussFilt(fftshift(real(ifftn(fftn(m1d).*conj(fftn(m2d))))),.1);
%         [~,shifts]=max3di(cc);
%         zsh=zsh+shifts(3)-ctr;
        fSh=FourierShift(n1,[0 0 zsh]);
        m1dsh=real(ifftn(fftn(m1d).*fSh));
        m1dr=rsRotateImage(m1dsh,r);
        err=sum((m1dr(:)-m2d(:)).^2);
        [P,t]=Simplex(err);
        zsh=P(1);
                r=P(2);

        if max(std(t.simp))<minAngle
            break;
        end;
        if mod(i,displayIters)==0
            subplot(2,2,1);
            imags(sum(m1dr,1));
            title(['ds ' num2str(ds1) '  iter ' num2str(i) '  shift ' num2str(zsh)]);
            subplot(2,2,2);
            imags(sum(m1dr,3))
            title(['rot ' num2str(r) '  err ' num2str(err)]);
            subplot(2,2,3);
            imags(sum(m1dr-m2d,1));
            subplot(2,2,4);
            imags(sum(m1dr-m2d,3));
            drawnow;
            disp([num2str([ds1 i]) '  ' num2str([zsh*ds1 r]) '  ' num2str(err)]);
%             disp(t.simp)
        end;
    end;
    zsh=zsh*ds1;
    rot=r;
end;