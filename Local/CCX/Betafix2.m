function ok=Betafix2(C)
% function ok=Betafix2(C)
% To get an empty contingency struct, call C0=Betafix();
% Fields of contingency C
% bPrice     Betafix price
% volume     Volume of combination

if nargin<1
    C=struct;

    C.bPrice.val=0; % all value fields may be n x 1 vectors.
    C.bPrice.limits=[0 1e4];
    C.bPrice.string='Betafix unit price';
    
    C.volume.val=0;
    C.volume.limits=[0 1e5];
    C.volume.string='Volume';
    
    ok=C; % return the struct
else
    p1=4000;  % price at vol=1
    p80=1500; % price at 80k
    k1=1/8e4;
    p81=2000; % price at 81k
    p500=1200;
    k2=1/5e5;
    
    % Make a function that decays from 1 at volume=1 to 0 at volume = 80k
    lowVol=C.volume.val<=8e4;
    lVols=C.volume.val(lowVol);
    f80=(exp(-k1*lVols)-exp(-k1*8e4)) ...
        /(exp(-k1)-exp(-k1*8e4));
    minPrice(lowVol)=p80+f80*(p1-p80);
    %     decays from 1 at volume =80k to 0 at volume = 500k
    hVols=C.volume.val(~lowVol);
    f500=(exp(-k2*hVols)-exp(-k2*5e5)) ...
        /exp(-k2*8e4)-exp(-k2*5e5);
    minPrice(~lowVol)=p500+f500*(p81-p500);
    
    ok=(C.bPrice.val>=minPrice') & (C.volume.val<=5e5); % not valid beyond 500k
    
end;