% deMakeJpegSequence
figure(1);
clf;
SetGrayscale;
fc=.1;

[sFile inDir]=uigetfile('*','Pick -s file');
cd(inDir);
label=[sFile(1:8) sFile(46:50)];

load(sFile);
cd ../JpegSequence
nim=size(s.imgs,3);
% get scaling
    img=(single(s.imgs(:,:,2))-s.dr1)./(s.br1-s.dr1);
    n=size(img);
    img=GaussFilt(img,fc);
    [img8 mul add]=imscale(Downsample(img,n/4),256,1e-4);
for i=1:nim
    img=(single(s.imgs(:,:,i))-s.dr1)./(s.br1-s.dr1);
        img=GaussFilt(img,fc);

    img=Downsample(img,n/4);
    img8=uint8(img*mul+add);
    imac(img8);
    outName=[label sprintf(' %02d',i) '.jpg'];
    title([outName num2str([mul add])]);
    drawnow;
    imwrite(rot90(img8,3),outName);
end;
hist(img(:)*mul+add,1000);