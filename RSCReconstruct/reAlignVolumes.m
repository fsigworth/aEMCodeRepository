function [vol2ali,tz,gamma,mirror]=reAlignVolumes(vols1,vols2,free)
% function vol2a=reAlignVolumes(vols1,vols2,testMirror)
% Do rotational and translational alignment of volumes, and return the
% aligned vols2.  The free parameter flags are:
% [shiftxy shiftz gamma mirror] with defaults [0 1 1 1]
% The rotation is tested just using the Z-projection of the volume.  There
% is code for doing a full 3D rotation about Z, but it's commented out.
% A soft mask is applied in line 29, msk=fuzzymask(n,3,0.35*n,.1*n);


% cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon96bo/mrc');
% vols1=ReadMRC('i22av01.mrc');
% vols2=ReadMRC('i22bv01.mrc');
% testMirror=1;
% vols2=MirrorX(vols1)+0.3*randn(size(vols2));
% figure(1);

if nargin<3
    free=[0 1 1 1];  % shiftxy shiftz gamma mirror
end;

n=size(vols1,1);
nv=size(vols1,4);
vol2ali=vols1;
msk=fuzzymask(n,3,0.35*n,.1*n);
gammaStep=[20];
gammaMin=2;
tz=zeros(nv,1);
tgamma=zeros(nv,1);
mirror=zeros(nv,1);

for iVol=1:nv
    v1=vols1(:,:,:,iVol);
    v1m=v1.*msk;
    v2=vols2(:,:,:,iVol);
    %     Get the z-projections
    p1=sum(v1m,3);
    %     p2=sum(v2,3);
    for iter=1:numel(gammaStep)
        %         Get the displacement.  Should be nonzero only in z
        cc=fftshift(real(ifftn(fftn(v1m).*conj(fftn(v2)))));
        [~,p]=max3di(cc,1);
        tz(iVol)=tz(iVol)+p(3);
        if ~free(1)
            p(1:2)=0;  % shift only in z
        end;
        %          p
        if free(2)
            v2sh=real(ifftn(fftn(v2).*FourierShift(size(v2),p)));
        else
            v2sh=v2;
        end;
        v2m=v2sh.*msk;
        %         Get the rotation from the projection
        p2=sum(v2m,3);
        ccs=zeros(2,1);
        gammas=zeros(2,1);
        for iMirror=1:1+(free(4)>0)
            %              iMirror
            if free(3)
                gamma=Simplex('init',0,gammaStep(iter));
                done=0;
                while ~done
                    p2r=rsRotateImage(p2,gamma);
                    cc=p2r(:)'*p1(:);
                    %                 v2r=rsRotateImage(v2m,gamma);
                    %                 cc=v2r(:)'*v1m(:);
                    %                 subplot(221); imags(p1);
                    % %                 subplot(222); imags(sum(v2r.*msk,3));
                    %                 subplot(222); imags(p2r);
                    %                 title(gamma);
                    %                 drawnow;
                    %             err=1e6*(gamma<-45 | gamma>45);
                    gamma=Simplex(-cc);
                    %                 disp([gamma cc/1e7])
                    done=Simplex('converged',gammaMin);
                end;
                gamma=Simplex('centroid');
                p2r=rsRotateImage(p2,gamma);
            else
                gamma=0;
                p2r=p2;
            end;
            ccs(iMirror)=p2r(:)'*p1(:);
            gammas(iMirror)=gamma;
            v2m=MirrorX(v2m);
            p2=MirrorX(p2);
            %             p2=MirrorX(p2);
        end;
        [~,iMirror]=max(ccs);
        %         iMirror
        gamma=gammas(iMirror);
        %         Rotate the volume
        if iMirror==2
            v2sh=MirrorX(v2sh);
        end;
        v2=rsRotateImage(v2sh,gamma);
    end;
    vol2ali(:,:,:,iVol)=v2;
    %     subplot(221); imags(sum(v1,3));
    %     subplot(222); imags(sum(v2,3));
    %     drawnow;
    mirror(iVol)=iMirror-1;
    tgamma(iVol)=gamma;
    
end;



