function [minPrice,ind]=Payer2(C)
% function [minPrice,ind]=Payer2(C)

if nargin<1  % return a dummy structure
%     C.isCancer.val=true;
    C.aPrice.val=0;
    C.bPrice.val=0;
    C.volume.val=0;
    C.aSupplyChain.val=true;
    C.aWarehouse.val=true;
    ok=C;  % return this struct
else % normal call
%     if C.isCancer
%  Test for prices up to limitPrice
    [minPrice,ind]=min(C.aPrice.val+C.bPrice.val); % find the lowest
end;
