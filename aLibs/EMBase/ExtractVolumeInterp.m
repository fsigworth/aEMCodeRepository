function out=ExtractVolumeInterp(m,cp,n);
% Get a smaller volume of size m from the 3D map m, with center point cp.
% It's like ExtractVolume except we do trilinear interpolation
% We don't do any clipping yet!

epsi=1e-4; % no interpolation necessary if fractional shift is smaller than this.
if numel(n)<3
    n=[1 1 1]*n(1);
end;
% m=map;
% n=32*[1 1 1];
% % m=zeros(200*[1 1 1],'single');
% cp=[119.52 130.44 140.33]/1.068;
% m(cpi(1),cpi(2),cpi(3))=1;
cp1=cp-floor(n/2)-1;
int=floor(cp1);

% Integer coordinates
ixs=int(1)+1:int(1)+n(1);
iys=int(2)+1:int(2)+n(2);
izs=int(3)+1:int(3)+n(3);

% Residual fractional coords
frac=cp1-int;
if all(abs(frac)<epsi)
    out=ExtractVolume(m,cp,n);
else
    %     frac=zeros(1,3);
    % frac=[-1 -2 -1];
    fx=frac(1);
    fy=frac(2);
    fz=frac(3);
    
    out=(1-fz)*( (1-fy)*( (1-fx)*m(ixs,iys,izs)  +fx*m(ixs+1,iys,izs)  )...
        +fy*( (1-fx)*m(ixs,iys+1,izs)+fx*m(ixs+1,iys+1,izs) ) )...
        +fz*( (1-fy)*( (1-fx)*m(ixs,iys,izs+1)+fx*m(ixs+1,iys,izs+1) )...
        +fy*( (1-fx)*m(ixs,iys+1,izs+1)+fx*m(ixs+1,iys+1,izs+1) ) );
    
    %  ShowSections(out);
end;
