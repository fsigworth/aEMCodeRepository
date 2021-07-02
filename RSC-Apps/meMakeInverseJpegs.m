% meMakeInverseJpegs.m
% Find all the relevant images in the given source directory. Apply the CTF
% inverse filter with the given compensation fraction, and also apply a
% Gaussian filter as desired.
% % miNames=f2FindInfoFiles;
if numel(miNames)<1
    disp(['No mi files found in ' pwd]);
end;
sourceDir='Merged_sm/';
sourcePattern='ms.mrc';
targetPattern='msi.jpg';
targetDir='JpegInv_3/';
unfiltTargetPattern='ms.jpg'; % if not empty, we write this too.
% unfiltTargetPattern=''; % if not empty, we write this too.

doWrite=1;

% compFraction=.2;
% disFc=1; % Gauss filter cutoff
% disHp=.01;
% clipThresh=1e-4;
% mulr=600;
% addr=160;

compFraction=.3;
disFc=.5; % Gauss filter cutoff
disHp=.014;
clipThresh=1e-4;
mulr=.7;
addr=160;


imags(zeros(256));
CheckAndMakeDir(targetDir,1);
nFound=0;
nErrs=0;
nt=numel(targetPattern);
for i=1:numel(miNames)
    mi=ReadMiFile(miNames{i});
    if ~isfield(mi,'padImageSize')
        mi.padImageSize=mi.imageSize;
    end;
    sourceName=[sourceDir mi.baseFilename sourcePattern];
    if ~exist(sourceName,'file')
        nErrs=nErrs+1;
        if nErrs<10
            disp([sourceName ' : file not found.']);
            continue;
        elseif nErrs==10
            disp('Continuing to search...');
            continue;
        else
            title(['searching ' num2str(i)]);
            drawnow;
            continue;
        end;
    end
    nErrs=0;
    ms=ReadEMFile(sourceName);
    outInvName=[targetDir mi.baseFilename targetPattern];
                msInv=rspCTFInverseFilter(ms,mi,compFraction);
%                 msInvDis=GaussFilt(msInv,disFc);
msInvFilt=GaussHP(GaussFilt(msInv,disFc),disHp);
%                 msInvDis=max(0,min(1,imscale(msInvFilt,1,clipThresh)));
                msInvDis=max(0,min(256,mulr*msInvFilt+addr));
% 
%                 nc=round(size(ms)*.7); % we'll use the center for computing clipping.
%                 mCtr=Crop(msInvDis,nc);
%                 npCtr=numel(mCtr);
% 
%                 [~,mulr,addr]=imscale(mCtr,1,clipThresh);
%                 msInvDis=max(0,min(1,msInvDis*mulr+addr));
% %                 
% %                 % Code copied from imscale
% %                 [mcSort,inds]=sort(mCtr(:));
% %         mnp=mCtr(round((npCtr-1)*clipThresh+1)); % lowest pixel value
% %         mxp=mCtr(round(npCtr-(npCtr-1)*clipThresh)); % highest pixel value
% %         
% %         msInvDisC=max(mnp,min(mxp,msInvDis));
        if numel(unfiltTargetPattern)>1
            outNfName=[targetDir mi.baseFilename unfiltTargetPattern];
            msNfDis=WriteJpeg(ms,outNfName,1e-3,doWrite); % take the default clipping
            imaga(msNfDis);
            pause(0.1);
        end;
        WriteJpeg(msInvDis,outInvName,0,doWrite); % no clipping
        imags(msInvDis);
        title([num2str(i),': ' outInvName],'interpreter','none');
        drawnow;
end;
