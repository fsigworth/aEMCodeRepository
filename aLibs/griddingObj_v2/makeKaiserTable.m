function w=makeKaiserTable(kernelsize,ov)
% Create the interpolated Kaiser table with oversampling factor ov
% the returned w array is ov x kernelsize in size.
% Modified from Fred's code - removed sinc mode

persistent w1;  % Cache for the table

nw=kernelsize;
mo=lower(mode(1));

    alpha=gridGetAlpha(kernelsize);
    ov=1024;  % oversampling of window for look-up

if (numel(w1)~=nw*ov) % if the cached table is absent or not the right size
    % Create the interpolation kernel
    epsi=1e-8;  % Crude way to avoid divide by zero in the sinc function.
    dw=1/ov;  % fractional increment.  We assume ov is even
    k=(-nw/2+dw/2:dw:nw/2-dw/2)';  % nw*ov space points
    w=gridLibKaiser(nw/2,alpha,k);  % Compute the Kaiser window
    w=w*ov/sum(w);  % normalize it.
    w1=zeros(ov,nw);
    % Make the 1D lookup table w1(i,:).  i=1 corresponds to a shift between -0.5 and
    % 0.5+dw, so we assign it the mean shift -0.5+dw/2.
    % i=ov corresponds to a mean shift of 0.5-dw/2.
    for i=1:ov
        w1(i,:)=w(ov-i+1:ov:nw*ov);
    end;
end;
w=w1'; % transpose to (nw x ov) size.
