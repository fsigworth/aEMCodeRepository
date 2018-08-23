function logi=BitsetToLogical(bits)
% Decode a bitset (packed bytes structure) into a logical array.
sz=bits.size;
b=bits.bytes;
logi=false(numel(b),8);
for i=1:8
    logi(:,i)=mod(b,2);
    b=idivide(b,2);
end;
logi(prod(sz)+1:end)=[];
logi=reshape(logi,sz);

