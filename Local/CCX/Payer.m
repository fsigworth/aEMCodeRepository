function ok=Payer(C)
% function ok=Payer(C)
% Check if contingencies are ok
%     ok=Payer(C);
% Return an empty contingencies struct (to derive names)
%     C0=Payer();

if nargin<1  % return a dummy structure
    C.isCancer=true;
    C.aPrice=0;
    C.bPrice=0;
    C.volume=0;
    C.aAltCompensation=0;
    C.aSupplyChain=true;
    C.aWarehouse=true;
    ok=C;  % return this struct
%     min 70,000; max 90,000
%     min vol 100, max 120

else % normal call
    if C.isCancer
        priceOk=C.aPrice+C.bPrice <= 5000 && C.volume <= 1.2e5;
    else
        priceOk=C.aPrice+C.bPrice <= 4000 && C.volume <= 1e6;
    end;
    
    compOk=C.aAltCompensation<200*C.volume;
    
    ok = priceOk && compOk && C.aSupplyChain && C.aWarehouse;
end;