imagePath='/Volumes/TetraData/EMWork/Scott/';
 imageName='11oct20b_20170729_04_4096x4096_bright_0.mrc';
imageName='carbon.mrc';
% imageName='floodbeam.mrc';
% imageName='10sep19a_a_00006gr_00022sq_v01_00002hl_v01_00003en.mrc';
inName=[imagePath imageName];
outName=[inName 'o'];
varFactor=2;
showSpectra=1;

% function ok=CompressMRCFile(inName, outName, varFactor, showSpectra) 
% % function CompressMRCFile(inName, outName, varFactor, imageMode);
% % Read an mrc file, compress the data and write out 
% if nargin<4
%     showSpectra=1;
% end;
% imageMode=1;
% if nargin<3
%     varFactor=8;  % about 3 bits/pixel
% end;

[m0 s0]=(ReadMRC(inName));
% Code for handling various datatypes; but we always write out float32.
% if isa(m0,'int16')
%     mrcMode=1;
% elseif isa(m0,'single')
     mrcMode=2;  % floating
% else
%     error('Unrecognized input file type');
% end;
m=double(m0);
n=size(m,1);

%%
f=double(Radius(n))/n;
sp2=CCDModelSpectrum(f,1);
sp2=sp2/(sp2(n/2+1,n/2+1));  % normalize to 1
fm=fftn(m);
fmw=fm./ifftshift(sqrt(sp2+.0));
mw=real(ifftn(fmw));  % pre-whitened
% Estimate the quantum
msk=fuzzymask(n,2,n*0.4,1)-fuzzymask(n,2,n*0.3,1);
% msk=1;
imgMean=mean(mw(:));
mskFm=fftshift(fmw).*msk;
imgVar=abs(mskFm(:)'*mskFm(:))/(msk(:)'*msk(:)*n^2);  % mean spectral density
quantum=imgVar/imgMean;

% Try compression
% sensitivity is 33 counts/e; assume quantal sens is 50 counts/e.
% varFactor=8;  % variance factor  % gives 3 pits/pixel, >100:1
% varFactor=4;  % 2 b/pixel, ~25:1
rsens=quantum;  % convert to counts
k=2*varFactor/sqrt(rsens);
msq=sign(mw).* floor(k*sqrt(abs(mw)));

figure(1);
SetGrayscale;
subplot(2,2,1);
imacs(BinImage(m,8));
title(inName,'interpreter','none');
axis off
subplot(2,2,2);
hist(msq(:)/k,1000);
title(['quantum = ' num2str(quantum)]);
xlabel(['sqrt(p)']);
drawnow;

% Do the reconstruction
% Find the center of the interval between the squares of
% msk/k (msk/k+1/k)
mss=msq/k;
mrq=sign(mss).*(abs(mss).^2+abs(mss)/k+0.5/(k^2));  % shift to center of bin

mr=real(ifftn(fftn(mrq).*ifftshift(sqrt(sp2))));

subplot(2,2,1);
imacs(BinImage(mr,8));
title(outName,'interpreter','none');
axis off;

subplot(2,2,3);
plot([sect(m) sect(mr)-sect(m)]);
legend('Input','Difference','location','east');
xlabel('x-position of central line');
drawnow;

% Write out the modified file.
fo=WriteMRCHeader(mr,s0.pixA,outName,1,-n/2*[1 1 0],mrcMode);
if mrcMode==2
    count=fwrite(fo,mr,'single');
elseif mrcMod==1
    mr=round(mr);
    count=fwrite(fo,mr,'int16');
end;
    fclose(fo);
ok=count>0;

% 
if showSpectra
    subplot(2,2,4);
spm=RadialPowerSpectrum(m);
spmr=RadialPowerSpectrum(single(mr));
freqs=0:1/n:1/2-1/n;
semilogy(freqs,[spm abs(spm-spmr)]);
legend('Input spectrum','Difference spectrum');
ylabel('Spectral density');
xlabel('Relative frequency');
drawnow;
end;

%%
if 1
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
            blockData=msq(x0:x1,y0:y1);  %% note max
            bd=blockData(:);
            x=min(bd):max(bd);
            h=hist(bd,x);
            [hs inds]=sort(bd,'descend');
            vals=x(inds);
            vals=vals(hs>0);
            hs=hs(hs>0);
            disp(dict);
            dict=huffmandict(vals,hs/sum(hs));
            disp(hcode);
            hcode=huffmanenco(bd,dict);
%             
%             bR(i,j)=bitRate(1);
%             disp([i j bitRate(1)]);
%             code=[code;uint8(codedBlock)];
        end;
    end;
    toc
    
end;

return
%%
% 
% n=4096;
% 
% mx=zeros(n,n);
% my=mx;
% 
% s1=2.5;
% s2=-.7;
% ct=n/2+1;
% 
% mx(ct,ct)=s1+s2+1;
% mx(ct+1,ct)=-s1;
% mx(ct+2,ct)=-s2;
% my(ct,ct)=s+1;
% my(ct,ct+1)=-s;
% 
% fx=fftn(fftshift(mx));
% % fy=fftn(fftshift(my));
% sx=sectr(fftshift(fx));
% % sz=sectr(fftshift(fx.*fy));
% % semilogy(abs([sx sz]));
% subplot(2,1,2);
% semilogy(abs(sx).^2.*sp0);

