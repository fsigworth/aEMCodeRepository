% k2CoincidenceLossSpectrum

% cd('/Volumes/Extreme1T/20181025/80Frames/mv')
% mvName='grid9_0006.tif';
% disp(['Reading ' mvName]);
% [mv,s]=ReadMovie(mvName);
% fprintf('%d frames\n',s.nz);
% fprintf('Removing outliers');
% mf=zeros(size(mv),'single');
% for i=1:s.nz
%     mf(:,:,i)=RemoveOutliers(mv(:,:,i));
%     fprintf('.');
% end;
% fprintf('\n');

% n=NextNiceNumber(max(s.nx,s.ny),5,8);
% disp(n);
% mp=zeros(n,n,s.nz,'single');
% means=zeros(s.nz,1);
% % Pad the frames
% w=SquareWindow([s.nx s.ny],32);
% for i=1:s.nz
%     imgf=mf(:,:,i);
%     means(i)=mean(imgf(:));
%     mp(:,:,i)=means(i)+Crop((imgf-means(i)).*w,n);
% end;
% m=sum(mp,3);
% disp('Spectra');
% tic;
% sps=RadialPowerSpectrum(mp,1);
% toc;
% spectrum=sum(sps,2);
% save spectrum sps spectrum
% 
% sp=spectrum;
% sp(1:4)=sp(5);

a=125.5;
b=.17;
d=2.2;

fs=(1:n/2)/n;

S=a*(1-b*(sinc(fs*d)).^2)';

plot([sp S]);

