function ok=AlphaLabs(C)
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
    C.aSupplyChainPenalty=0;
    C.aWarehouse=1;
    ok=C;
else   
    
    
%     20,000 per dose if vol < 100
%     reduce by 20% with warehouse
%     >= betafix price
%     supply chain penalty 300k
%     Total revenue has to be at least 2M
    totalRevenueOk=C.aPrice*C.volume+C.aAltCompensation>2e6;
    
    % Thresholds for price
p1=2e4;
    % Lower price if there's a warehouse in Zurich
    if C.aWarehouse
        p1=p1*0.8;
    end;
    
    ok= totalRevenueOk ...
        && ~C.aSupplyChainPenalty ...
        && C.aPrice>=p1 ...
        && C.aPrice >= 1 * C.bPrice;
    
end;
