% TestCompNoise
% Pick a file pair and compute the compensation noise.
% See how well it eliminates the acf peak.
% We test two functions: meMakeFixedPatternCompNoise and ...Auto
%  It seems that the Auto function works about as well.
CCDCPE=16;

filter='*.mrc;*.dm3';
[fname fpath]=uigetfile(filter,...
    'Select three micrograph files for noise comp','multiselect','on');
cd(fpath);
if isa(fname,'numeric') % user clicked cancel
    return;  % exit from the program.
elseif isa(fname,'char')  % simple string
    error('Need 3 files');
else
    numfiles=numel(fname);  % cell array
end;
%%
disp('Reading images');
[m pixA doses]=meReadImages(fname,CCDCPE);
m=mePreWhiten(m);
disp('Making comp noise from images 1 and 2');
% [CompNoise cc]=meMakeCompNoise(m(:,:,1),m(:,:,2),prod(doses(1:2)));
[CompNoise cc]=meMakeCompNoiseAuto(m(:,:,1),m(:,:,2));

% disp('Storing the noise');
% save CompNoise CompNoise

disp('testing')

% pa=fileparts(which('CCDMakeFixedPatternRefCCF'));
% load([pa '/FixedPatternCCF.mat']);  %loads ccf0
% % nc=size(ccf0,1);
%
% ccq=prod(doses)*fftshift(real(ifftn(abs(fftn(CompNoise)).^2)));
%
figure(1); SetGrayscale;
subplot(3,2,2);
imacs(Crop(cc,nc));
% % plot(sect(Crop(ccq,nc))*prod(doses));
title('noise acf');

% subplot(2,2,2);
% % imacs(m(:,:,2));
% plot(sect(ccf0));
% title('ccf0');

subplot(3,2,1);
imacs(Crop(cc,nc));
plot(sect(Crop(cc,nc)));
title('Raw CC');
drawnow;
for i=1:2
    cc2=fftshift(real(ifftn(fftn(m(:,:,i)+doses(i)/doses(1)*CompNoise).*...
        conj(fftn(m(:,:,i+1)-doses(i+1)/doses(2)*CompNoise)))));
    subplot(3,2,2*(i+1));
    % imacs(Crop(cc,nc)-2*ccf0);
    imacs(Crop(cc2,nc));
    subplot(3,2,2*(i+1)-1);
    plot(sect(Crop(cc2,nc)));
    % plot(sect(Crop(cc2,nc)));
    title('Comp CC');
    drawnow;
end;

