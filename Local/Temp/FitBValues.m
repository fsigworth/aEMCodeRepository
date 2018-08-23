% FitBValues.m
% Estimate B values in K2 images
%
fileStart=8;
fileEnd=8;

nq=1;  % nq x nq quadrants
nxq=floor(3840/nq);

Pa=struct;  % fitting parameters
Pa.lambda=EWavelength(200); % lambda in angstroms
Pa.Cs=2;
 Pa.defocus=1:.2:3;
% Pa.defocus=7:.2:9;  % for image 2
Pa.deltadef=-.1:.05:.1;  % astigmatism
Pa.theta=0:.3:pi/2;
Pa.alpha=.02;  % This can also be fitted, for use with phase plate.
Pa.B=20:10:60;  % This is fitted if a vector of values is given.

nu=1024;
ds=1;
fRange=[.05 .3];


figure(1);
SetGrayscale;
xlabel('Intensity');
ylabel('B');


cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1');
d1=dir('Info/');
Ps=struct;
Intens=zeros(fileEnd,nq,nq);
Mis={};
for fi=fileStart:fileEnd
    %
    dnm=d1(fi+12).name;
    mi1=ReadMiFile(['Info/' dnm]);
    Mis{fi}=mi1;
    disp(dnm);
    imgName=mi1.imageFilenames{1};   %  Image 2
    m10=ReadEMFile([mi1.imagePath imgName]);
    
    if fi==fileStart   % first image, initialize variables
        opts=struct;
        opts.blockSize=nu/ds;
        fn=RadiusNorm(opts.blockSize)/ds; % normalized frequency
        fw=fn/mi1.pixA;  % actual frequency
        
        % Parameters for carbon model
        a0=50;
        k4=.034;
        a4=2080;
        expo=3.5;
        
        mtf2=CCDModelMTF(fn,5).^2;
        opts.carbonModel=mtf2.*((k4^expo)./(fw.^expo+k4.^expo)+a0/a4);
    end;
    
    for qi=1:nq^2 % quadrant index
        qix=mod(qi-1,nq)+1;
        qiy=floor((qi-1)/nq)+1;
        n=nxq;
        opts.title=[dnm(1:8) '  ' num2str([qix qiy])];
        m1=m10((qix-1)*nxq+1:qix*nxq,(qiy-1)*nxq+1:qiy*nxq);
        
        m2=Downsample(m1,n/ds);
        Intens(fi,qix,qiy)=mean(m2(:));
        
        % ---------fitting is done here ------------
        [P, cc]=FitCTFk2(m2,Pa,mi1.pixA,ds,fRange,opts);
        
        disp([dnm num2str([qix qiy])]);
        P
        disp(' ');
        disp(' ');
        if fi==fileStart && qi==1
            Ps=P;
        end;
        Ps(fi,qix,qiy)=P;
    end;
end;

return
%%
% save Ps1Quadrants Ps Intens
