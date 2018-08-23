function ok=Carpathia(C)
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
    C.aSupplyChainPenalty=false;
    C.aWarehouse=true;
    ok=C;  % return this struct
%     min 70,000; max 90,000
%     min vol 100, max 120

else % normal call
        sumPrice=C.aPrice+C.bPrice;
        totalPrice=sumPrice*1.01;
        priceOk=totalPrice <= 9e4 && totalPrice>= 7e4;
        volOk=C.volume>=100 && C.volume <=120;
    ok = priceOk && volOk && ~C.aSupplyChainPenalty && C.aWarehouse;
end;