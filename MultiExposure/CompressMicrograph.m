% CompressMicrograph.m

fname=['/Volumes/TetraData/EMWork/Hideki/110628/HS12slot4_GluR2-GFP/' ...
        '001u2000.dm3'];
countsPerE=1;
oversampling=1;
scale=4*oversampling/countsPerE;
    m=ReadEMFile(fname);
m=Crop(m,1024);  % sub-problem
    n=size(m,1);

me=mean(double(m(:)));
mw=mePreWhiten(single(m))+me;
%%

msq=floor(sqrt(mw*scale));
h=hist(msq(:),200);
% plot(xs,(h));

%%
nb=256;  % block size
blocks=ceil(n/nb);

bR=zeros(blocks,blocks);
disp('Huffman encoding');
code=[];
tic
for j=1:blocks
    y0=(j-1)*nb+1;
    y1=min(j*nb,n);
    for i=1:blocks
        x0=(i-1)*nb+1;
        x1=min(i*nb,n);
        blockData=msq(x0:x1,y0:y1);
        [codedBlock bitRate]=huff03(blockData(:),0);
        bR(i,j)=bitRate(1);
        disp([i j bitRate(1)]);
        code=[code;uint8(codedBlock)];
    end;
end;
toc
return


msq1=round(max(msq(1:1e5),0)');
[y br]=huff03(msq1,0);
whos y
br
yi=uint8(y);
dm=huff03(double(yi));
