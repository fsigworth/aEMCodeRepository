function [ctfs spectra]=meFitCTFs(m,pixA,defoci,doRemoveSpots,options)
% function [ctfs spectra]=meFitCTFs(m,pixA,defoci,doRemoveSpots)
% Given the array of images m, create an array of ctf parameter structures.
% Also return the averaged spectra in the cell array spectra:
% spectra{i}.df        % frequency step in A^-1
% spectra{i}.spectrum  % square 2D spectrum, f=0 in center
% defoci is a vector of the approximate defocus values.
% options have the fields
% spFromWholeMicrograph: use entire micrograph for spectrum estimate.
% lowB: use lower B values (half of normal)
% kV: use the given kV value instead of 200

if nargin < 4
    doRemoveSpots=0;
end;
if nargin < 5
    options=struct;  % empty struct
end;

validOptions={'kV','lowB','spFromWholeMicrograph'};
optionDefaults=[200 0 0];
options=ScanOptionsStruct(options,validOptions,optionDefaults);

ctFitPars=meSetCTFitPars(defoci,pixA,options.kV);

if isfield(options,'lowB') && options.lowB
    for i=1:numel(ctFitPars)
        ctFitPars(i).B=ctFitPars(i).B/2;
    end;
end;

% Do the CTF fitting
[nx ny nim]=size(m);
n=[nx ny];
for i=1:nim
    if defoci(i)>=1.5 % large defocus; fit CTF with binned image
        ds=2;
        maxres=10;
        minres=100;
%         minres=200;
    else            % small defocus: fit CTF with full image
        ds=1;
        maxres=max(pixA*3, 5);  % fit out to max 2/3 of Nyquist
        minres=50;
    end;
    ns=n/ds;
    [mc fmc]=Downsample(m(:,:,i),ns);  % reduce image size by ds
%     
%     figure(4); clf;
%     SetGrayscale;
%     imacs(mc);
%     drawnow;

    if doRemoveSpots
%         ns1=n/ds1;
        nx=max(ns);
        if nx>min(ns)  % rectangular image
            win=SquareWindow(ns);
            mcpad=Crop(win.*mc,nx);  % downsampled, padded image
            fmc1=fftn(mcpad);
        else
            fmc1=fmc;
        end;
        fmcs=fftshift(fmc1);  % zero-centered
        ds1=max(ds,2);
        nr=nx/ds1;
        if ds1>ds
            fmcs=Crop(fmcs,nr);
        end;
        minr=nr*pixA*ds1/70;
        threshSD=5;
        display=1;
%         disp('RemoveSpots');
        [spc pmask]=RemoveSpots(abs(fmcs).^2,minr,threshSD,display);
        % pmask is now zeros at the spots, and square.
        xpmask=Crop(pmask,nx,0,1);  % fill the outside with ones
        rmask=Downsample(xpmask,ns);  % Convert to rectangle
        mc=real(ifftn(ifftshift(rmask).*fmc));
    end;
    
%     disp(['CTF fitting ' num2str(i)]);
    [P c sp]=FitCTF(mc,ctFitPars(i),pixA*ds,maxres,minres,1,options);
    drawnow;
    spectra{i}.spectrum=sp;
    spectra{i}.df=1/(size(sp,1)*pixA*ds);
    P.res=pixA;  % Change it back to the original pixel size
    if i==1
        ctfs=P;
    else
        ctfs(i)=P;
    end;
end;


