function [m, pixA, doses]=meReadImagesNorm(fnames, cpe, flipFirst,overridePixA, weights, skipOutliers)
% function [m pixA doses]=meReadImagesNorm(mi, cpe, flipFirst,overridePixA, weights)
% function [m pixA doses]=meReadImagesNorm(fnames, cpe, flipFirst,overridePixA, weights, skipOutliers)
%
% Read in the stack of images m (n x n x nims)
% This is the normalized version, where the image values are the fractional
% contrast.
% fnames is a cell array of file names, or else it is the mi struct, from
% which the filenames are constructed.  cpe is counts per electron, which
% is read from mi if present; otherwise the default is 16.  pixA is
% angstroms per pixel read from the last image, unless overridePixA > 0.
% If pixA can't be read from the micrograph, the mi.pixA value is used.
% doses is an array of doses in electrons per square angstrom, computed
% from masked images, cpe and pixA. The images are
% returned scaled as fractional variation from the doses value and
% have zero mean.
% Images are padded to NextNiceNumber(n,5,4) to handle K2 sizes.
% The optional weights argument is an array 1 x nfiles.  Any zero entry
% means the corresponding file will be skipped.
% The returned images m have had outliers removed and are windowed at the
% edges (square window, n/windowFraction pixels) and normalized to unity contrast.
beamEdgeA=300; % amount by which the threshold mask is extended.

if isa(fnames,'char') % A simple char array; convert to cell
    q=cell(1);
    q{1}=fnames;
    fnames=q;
end;
if nargin<2 || cpe==0
    cpe=16;  % counts per electron
end;

if isa(fnames,'struct')  % we passed the mi structure
    mi=fnames;
    nim=numel(mi.imageFilenames);
    fnames=cell(nim,1);
    for i=1:nim
        fnames{i}=[mi.basePath mi.imagePath mi.imageFilenames{i}];
    end;
    if isfield(mi,'cpe')
        cpe=mi.cpe;
    end;
end;

nfiles=numel(fnames);

if nargin<3   % Check if we have to flip the first (aligned) image
    flipFirst=0;
end;
if nargin<4
    overridePixA=0;
end;
if nargin<5 || all(weights==0)
    weights=ones(1,nfiles);
end;

if nargin<6
    skipOutliers=0;
end;
windowFraction=64;

% allocate arrays
doses=zeros(1,nfiles);
win=[];

% Load all the images and compute the doses.
for findex=1:nfiles
    if weights(findex)
        name=fnames{findex};
        if numel(name)>0
            name=CheckForImageOrZTiff(name);  % possibly pick up a compressed file.
            disp(name)
            [m0, pixA]=ReadEMFile(name);
            if overridePixA
                pixA=overridePixA;
            end;
            if pixA<=0.1
                pixA=mi.pixA;
            end;
            if findex==1
                switch flipFirst
                    case -1
                        m0=m0';
                    case {1,3}
                        m0=rot90(m0,flipFirst);
                    case 0
                    otherwise
                        error('unrecognized value for flipFirst');
                end;
            end;
            if skipOutliers
                m0=double(m0);
            else
                m0=double(RemoveOutliers(m0));
            end;
            n0=size(m0);
            n1=NextNiceNumber(n0,5,8);
            if any(n1>n0) % we need to resize
                me=mean(m0(:));
                m0=Crop(m0,n1,0,me);
                n0=n1;
            end;                
            if numel(win)==0  % need to create the window
                win=SquareWindow(n0,n0(1)/windowFraction);
                m=single(zeros([n0 nfiles])); % output array
            end;
            %         if pixA==0
            %             error('pixA value is zero.');
            %         end;
            % Window and remove mean. Convert the final images to single to save space.
            if skipOutliers
                win2=win;
            else
                mf=GaussFilt(m0,.1);
                win1=mf>median(mf(:))/2;
                if sum(win1)<.999*numel(win1)  % significant masking to be done
                    msk=GaussFiltDCT(win1,0.2*pixA/beamEdgeA)>.95;
                    win2=win.*msk;
                else
                    win2=win;
                end;
            end;
            m1=m0.*win2/(cpe*pixA^2);  % m1 is in units of electrons
            dose=sum(m1(:))/sum(win2(:)); % mean value of windowed image
            % remove the mean and store as e / A^2
            m(:,:,findex)=single(m1/dose-win2);  % Normalized, zero-mean image
            doses(findex)=dose;
        end;
    end;
end;
