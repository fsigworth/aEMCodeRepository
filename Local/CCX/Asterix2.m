function ok=Asterix2(C)
% function ok=Asterix2(C)
% Vectorized function, allows vectors for C.field.val variables.
% Note that each val field has to be the same size vector (including
% logical variables).
%
%  To get an initialized struct use this call: C0=Asterix2();

if nargin<1 % return the initialized structure
    C=struct;

    C.aPrice.val=0; % all value fields may be n x 1 vectors.
    C.aPrice.limits=[0 1e4];
    C.aPrice.string='Asterix unit price';
    
    C.bPrice.val=0;
    C.bPrice.limits=[0 1e4];
    C.bPrice.string='Betafix unit price';
    
    C.volume.val=0;
    C.volume.limits=[0 1e5];
    C.volume.string='Volume';
    
    C.aSupplyChain.val=true;
    C.aSupplyChain.limits=[true true];
    C.aSupplyChain.string='Asterix supply chain is secure';
    
    C.aWarehouse.val=true;
    C.aWarehouse.limits=[true true];
    C.aWarehouse.string='Payer provides Zurich warehouse for Asterix';
    
    ok=C;  % returned variaable is the struct
else   
    % Total revenue has to be at least 10M
    totalRevenueOk=(C.aPrice.val .* C.volume.val) > 5e6;
    
    % Thresholds for price
    lowVolume=C.volume.val<1e5;
    p1(lowVolume)=2000;
    p1(~lowVolume)=500;

    % Lower price if there's a warehouse in Zurich
    p1(C.aWarehouse.val)=p1(C.aWarehouse.val)*0.8;
    
    ok= totalRevenueOk ...
        & C.aSupplyChain.val ...
        & C.aPrice.val>=p1' ...
        & C.aPrice.val >= 0.8 * C.bPrice.val;
    
end;
