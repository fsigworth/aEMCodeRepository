function ok=Asterix(C)
% function ok=Asterix(C)
%  to get empty struct: C0=Asterix();
%
% We require values from these fields of C
% aPrice  Asterix price
% bPrice  Betafix price
% volume  Volume of combination
% aAltCompensation  Asterix alternative compensation
% aSupplyChain  (boolean) Supply chain is secure
% aWarehouse    (boolean) Payer provides Zurich warehouse

if nargin<1
    C.aPrice=0;
    C.bPrice=0;
    C.volume=0;
    C.aAltCompensation=0;
    C.aSupplyChain=1;
    C.aWarehouse=1;
    ok=C;
else   
    % Total revenue has to be at least 10M
    totalRevenueOk=C.aPrice*C.volume+C.aAltCompensation>1e5;
    
    % Thresholds for price
    if C.volume<1e3
        p1=200;
    else
        p1=50;
    end;
    
    % Lower price if there's a warehouse in Zurich
    if C.aWarehouse
        p1=p1*0.8;
    end;
    
    ok= totalRevenueOk ...
        && C.aSupplyChain ...
        && C.aPrice>=p1 ...
        && C.aPrice >= 1.1 * C.bPrice;
    
end;
