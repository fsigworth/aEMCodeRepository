% meProcessFilmImage.m
interactive=1;
dose=20;
defocus=3;
mcDS=2;
FitCTF=1;
writeInfo=1;
imagePath='/Volumes/TetraData/EMWork/ip3r/';
% imagePath=[imagePath 'MicrographsClosed/'];
imagePath=[imagePath 'MicrographsOpen/'];
cd(imagePath);

if interactive
    d=[];
    [files uiPath]=uigetfile('*','Select an image','');
    imagePath=uiPath;
    if isnumeric(files);
        return
    end;
    if ischar(files)
        nim=1;
        d.name=files;
    else
        mi.imageFilenames=files;
        nim=numel(files);
        for i=1:nim
            d(i).name=files{i};
        end;        
    end;
   
else  % process all files in the directory
    d=dir(imagePath);
    nim=numel(d);
end;

procPath=imagePath;
infoPath=procPath;


for ind=1:nim
    names=d(ind).name;
    if numel(names)>3 && strcmp(names(1:3),'dwc') % it's an image
        files=d(ind).name;
        
        
        mi=meCreateMicrographInfoStruct10;
        mi.imageFilenames{1}=files;
        
        
        [pa1 nm ex]=fileparts(mi.imageFilenames{1});
        imNumber=sscanf(nm,'dwc%d');
        
        if imNumber<1700  % closed state
            pixA=2.258;
            disp('Closed-state image');
        else
            pixA=2.8;
            disp('open-state image');
        end;
        
        
        mi.baseFilename=nm;
        disp(char(files));
        
        
        infos={};
        infos{1}=mi;
        startImage=1;
        
        nmi=numel(infos);
        %%
        i=1;
        disp(['>>Image index ' num2str(i)]);
        mi=infos{i};  % pick up the structure
        mi.imagePath=imagePath;
        mi.procPath=procPath;
        mi.infoPath=infoPath;
        
        % Construct the complete image file names
        nim=numel(mi.imageFilenames);
        fnames=cell(1,nim);
        figure(1);
        disp('Processing image:');
        fname=[mi.imagePath mi.imageFilenames{1}];
        %         disp(fnames{j});
        disp(fname);
        %%
        % Read the images
        %         [m mi.pixA mi.doses]=meReadImages(fnames);
        m1=single(ReadEMFile(fname));
        [nx ny]=size(m1);
        m1=m1(nx-ny+1:nx,:);
        n0=size(m1);
        mi.imageSize=n0;
        % get rid of saturated values
        m1(m1>250)=0;
        m1(m1<7)=0;
        nz=sum(m1(:)==0);
        %%
        me=sum(m1(:))/sum(m1(:)>0);
        m1(m1==0)=me;
        m1=m1-me;
        imacs(m1);
        
        
        
        % Window and remove mean. Convert the final images to single to save space.
        mi.pixA=pixA;
        mi.doses=dose;
        % remove the mean and store as e / A^2
        %         m(:,:,findex)=single(m1-win*dose);
        %         doses(findex)=dose;
        
        
        sz=size(m1);
        mi.nPixels=sz(1:2);
        
        %%
        %         operate with a pre-whitening filter
        h=exp(RadiusNorm(mi.imageSize)*.4);
        m=real(ifftn(fftn(m1).*fftshift(h)));
        %%
        % Fit the ctfs
        %         ncts=numel(mi.ctf);
        ncts=1;
        ctfOptions=[1 1];
        [mi.ctf mi.spectra]=meFitCTFs(m,mi.pixA,defocus,0,ctfOptions);
        defocus=mi.ctf.defocus
        %%
        % Align the images
        mi.mergeMatrix=eye(3);
        %%
        % Store the CTF and alignment information
        if writeInfo
            save([infoPath mi.baseFilename 'mi.mat'],'mi');
        end;
        
        %%
        % Create the composite image mc
        [effctf mc mts]=meCombineImages(m,mi.pixA,mi.ctf,mi.mergeMatrix,mi.doses,mcDS,0);
        mc=-mc;  % go back to reversed contrast
        %%
        LFAmp=.4;
        [mf c T]=meFilterImageFlatLF(mc,mi,LFAmp);
        figure(3);
        plot(sectr(c.*T));
        figure(2); SetGrayscale;
        mfg=GaussFilt(mf,.15);
        mfg1=GaussHP(mfg,.002);
        s=std(mfg1(:));
        mfgs=uint8(128+mfg1*256/s*.2);
        imac(mfgs);
        title([mi.baseFilename 'lfgs.jpg']);
        
        ofile1=[mi.procPath mi.baseFilename 'lf.jpg'];
        imwrite(uint8(imscale(mf)),ofile1);
        ofile2=[mi.procPath mi.baseFilename 'lfgs.jpg'];
        disp(['Writing ' ofile1]);
        imwrite(mfgs,ofile2);
        
    end; % if dwc file
end;  % for


