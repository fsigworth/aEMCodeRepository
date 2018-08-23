% MakeFixedPatternRef.m

% Pick up beam-noise images to make the ccf
path='/Users/fred/EMWork/Liguo/liguo-ccd-summer2010/';
fnames=[{'noise1.dm3'},{'noise2.dm3'}];
n=4096;
m=single(zeros(n,n,2));
for findex=1:2
    
    name=char(fnames(findex));
    name=[path name];
    disp(['Reading file: ' name]);
    [m1 pixA]=ReadEMFile(name);
    m(:,:,findex)=RemoveOutliers(m1);

    end;
    
%%
m=mePreWhiten(m);

ccf0=CCDMakeFixedPatternRefCCF(m);
figure(1); SetGrayscale;
subplot(2,2,1); imacs(ccf0);
subplot(2,2,2); plot(sect(ccf0));
