% PreWhitenParticles
% Read boxed "non-particles" and create a noise model.  Use this to filter
% micrographs.

writeOut=0;
noiseModelFcn='NoiseModel1';

mode=1;  % Operate on merged images
%  mode=3; % 110724 data
switch mode
    case 1
        pa='/Volumes/TetraData/EMWork/Hideki/110628/';
        boxPostName='mfBlank.box';
        imgPostName='m.mrc';
        imgPath='Merged/';
        boxPath='box_blanks/';
        preNames={'013';'014';'015';'039'};
        outputPreName='MergedWhitened/';
        miPath='mi/';
    case 3
        pa='/Volumes/TetraData/EMWork/Hideki/110724/';
        boxPath='box_Blank/';
        imgPath='Merge/';
        miPath='mi/';
        
        boxPostName='Blank.box';
        imgPostName='m.mrc';
        preNames={'19','39','59','79','119'};
        outputPreName='MergedWhitened/';
    otherwise
        error('Unrecognized mode');
        
end;

% Pick up stacks of non-particles
nstacks=numel(preNames);
doPause=0;
nb=64;
imgscale=1;
display=1;
sps=zeros(nb/2,nstacks);
figure(1);
SetGrayscale;

for j=1:nstacks;
    preName=preNames{j};
    imgname=[pa imgPath preName imgPostName];
    boxname=[pa boxPath preName boxPostName];
    [stack,pts,pixA]=ExtractBoxes(imgname, boxname, nb, imgscale, 1);
    title([preName imgPostName]);
    load([pa miPath preName 'mi.mat']);
    if j==1
        mis=mi;
    else
        mis(j)=mi;
    end;
    sp=RadialPowerSpectrum(stack);
    sps(:,j)=sp;
    if doPause
        pause;
    else
        drawnow;
    end;
end;
disp('loaded')

%%
figure(2);
co=get(gca,'colororder');

niters=800;

df=1/(pixA*nb);
freqs=(0:df:(nb/2-1)*df)';  % frequencies in spectrum of a single box
freqsx=(0:df/10:(nb/2-1)*df)';  % 10x oversampled frequencies


% Get the effective CTFs
effctf=[];
for j=1:nstacks
    [c effctf(:,j)]=meComputeMergeCoeffs(freqs, mis(j).ctf, mis(j).doses);
    [c effctfx(:,j)]=meComputeMergeCoeffs(freqsx, mis(j).ctf, mis(j).doses);
end;

% fit the noise model

% Initial parameters:
q=mean(sps(nb/4,:));

%  af1 af2 ag    sigma bf s0    f1exp f2exp
p=[ 5*q 5*q 25*q .007 100   q     1.5   2 ];
ac=[ 1   1    1    1   0    1     1    1  ];
p=Simplex('init',p,ac);
for i=1:niters
    [spec shot]=eval([noiseModelFcn '(freqs,p)']);
    model=repmat(spec,1,nstacks).*effctf.^2+...
        repmat(shot,1,nstacks);
    
    d=(model(:)-sps(:));
    err=d'*d;
    
    p=Simplex(err,p);
    if mod(i,10)==1
        set(gca,'colororder',co);
        plot(freqs,sps,'.');
        hold on
        plot(freqs,model,'-');
        hold off;
        legend(cellstr(preNames));
        title(i);
        xlabel('Spatial frequency, A^{-1}');
        drawnow;
    end;
    
end;
p=Simplex('centroid');

ok=input('Ready to filter the images? ','s');
if lower(ok)=='y'
    
    %%  Having made the model, do the filtering
    n=2048;
    f2d=Radius(n)/(n*pixA);
    
    % [coeffs effct]=meComputeMergeCoeffs(f2d,mi.ctf,mi.doses);
    % T=(ag*exp(-f2d.^2/(2*sigma^2)).*effct+s0)/s0;
    % Ti=1./sqrt(T);
    % save Ti Ti;
    
    disp('Filtering...');
    figure(1);
%     ctr=n/2+1;
%     H=NoiseModel1(f2d,p);
%     
%     gaussh=exp(-f2d.^2/(2*sigma^2));
%     f1h=(.01./f2d).^f1exp;
%     f1h(ctr,ctr)=f1h(ctr+1,ctr);
%     f2h=(.01./f2d).^f2exp;
%     f2h(ctr,ctr)=f2h(ctr+1,ctr);
%     H=(ag*gaussh+af1).*f1h+af2*f2h;  % Model that is multiplied by CTF
    
    d=dir([pa imgPath]);
    for i=1:numel(d)
        name=d(i).name;
        q=regexp(name,'m.mrc');
        if numel(q)>0  % contains 'm.mrc' in the name
            basename=name(1:q-1);
            [m pA]=ReadEMFile([pa imgPath name]);
            m=m-mean(m(:));
            n=size(m,1);
            load([pa miPath basename 'mi.mat']);
            if mi.version==9
                n0=mi.nPix(1);
            else
                n0=mi.imageSize(1);
            end;
              % Write out the noise model and parameters
            mi=meStoreNoiseModel(p,noiseModelFcn,mi);
            save([pa miPath basename 'mi.mat'],'mi');

            Ti=meGetNoiseWhiteningFilter(mi,n);
%             
%             pixA=mi.pixA*n0/n;
%             f2d=Radius(n)/(n*pixA);
%             
%             % Compute the effective ctf and inverse filter
%             [coeffs effct]=meComputeMergeCoeffs(f2d,mi.ctf,mi.doses);
%             [spec shot]=meEvalNoiseModel(f2d,mi);
% %             [spec shot]=eval([noiseModelFcn '(f2d,p)']);
%             T=(spec.*effct.^2+shot)./shot;
%             Ti=1./sqrt(T);
            subplot(1,2,1);                         %%
            plot(sectr(f2d),sectr(Ti));             %%
            axis([0 inf 0 1.02*max(sectr(Ti))]);    %%
            title('Prewhitening filter');           %%
            xlabel('Spatial frequency, Å^{-1}');    %%
            mf=real(ifftn(fftn(m).*ifftshift(Ti)));
            subplot(1,2,2);                         %%
            imacs(BinImage(mf,4));
            axis off;
            disp(name);
            title(name);
            drawnow;
            if writeOut
                WriteMRC(mf,pA,[pa outputPreName basename 'mw.mrc']);
                imwrite(uint8(imscale(mf)),[pa outputPreName basename 'mw.jpg']);
            end;
            %         basename
            %         pause
        end;
    end;
    disp('done.');
end;