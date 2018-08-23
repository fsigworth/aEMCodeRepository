function bits=BitsetFromLogical(logi)
% Pack an array of booleans into a bitset (an array of bytes), saving space.
bits.size=size(logi);
n=numel(logi);
nBytes=ceil(n/8);
logi=logi(:);
logi(end+1:nBytes*8)=false;  % Extend to a multiple of 8 elements
logi=reshape(logi,nBytes,8);
b=single(logi(:,1));
factor=1;
for i=2:8
    factor=factor*2;
    b=b+single(logi(:,i))*factor;
end;
bits.bytes=uint8(b);
