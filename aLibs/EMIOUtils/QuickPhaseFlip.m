% QuickPhaseFlip

    % Have the user select some info files
    [names, pathName]=uigetfile('*.*','Select image files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(names)
        names={names};
    end;
    cd(pathName);
%%



pa=struct;
pa.lambda=EWavelength(200);
pa.defocus=4:1:10;
pa.deltadef=-2:.4:2;
pa.theta=[-90:30:90]*pi/180;
pa.alpha=.02;
pa.B=100;
pa.Cs=2;



    
    for i=1:numel(names)
        [m,pixA,ok]=ReadEMFile(names{i});
        pixA=1.3;
        if ok
            m=sum(m,3);  % sum any movie
            n=size(m);
            figure(1);
            P=FitCTF(m,pa,pixA,10,50);
            P
            %%
            c=-ifftshift(CTF(size(m),P));
            k=.4;
            bInv=50;
            fs=ifftshift(RadiusNorm(n));
            boost=(.02+fs)./((fs/.06).^2+1);
            fm=BinImage(real(ifftn(fftn(m).*c.*boost./(k+c.^2))),4);
%             fm=SharpFilt(fm,.4);
%
figure(2);
SetGrayscale;
% imags(GaussFilt(m,.1));
imac(uint8(imscale(fm,256,.003)));
WriteJpeg(fm,[names{i}(1:end-4) '.jpg']);
        end;
    end;
    