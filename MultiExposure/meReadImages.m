function [m pixA doses]=meReadImages(fnames, cpe, flipFirst,overridePixA, weights)
% function [m pixA doses]=meReadImages(fnames, cpe, flipFirst,overridePixA, weights)
% function [m pixA doses]=meReadImages(mi, cpe, flipFirst,overridePixA, weights)
% 
% Read in the stack of images m (n x n x nims)
% fnames is a cell array of file names, or else it is the mi struct, from
% which the filenames are constructed.  cpe is counts per electron
% (default is 16).  pixA is angstroms per pixel for last image.  doses is
% an array of doses in electrons per square angstrom. The images are
% returned scaled in units of e/Å^2 and have zero mean.
% The optional weights argument is an array 1 x nfiles.  Any zero entry
% means the corresponding file will be skipped.

if isa(fnames,'char') % A simple char array; convert to cell
    q=cell(1);
    q{1}=fnames;
    fnames=q;
elseif isa(fnames,'struct')  % we passed the mi structure
    mi=fnames;
    nim=numel(mi.imageFilenames);
    fnames={};
    for i=1:nim
        fnames{i}=[mi.basePath mi.imagePath mi.imageFilenames{i}];
    end;
end;

nfiles=numel(fnames);

if nargin<2
    cpe=16;  % counts per electron
end;
if nargin<3   % Check if we have to flip the first (aligned) image
    flipFirst=0;
end;
if nargin<4
    overridePixA=0;
end;
if nargin<5
    weights=ones(1,nfiles);
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
            disp(name)
            [m0 pixA0]=ReadEMFile(name);
            if overridePixA
                pixA0=overridePixA;
            end;
            pixA=pixA0;
            if pixA0<=0
                pixA0=1;
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
            m0=double(RemoveOutliers(m0));
            n0=size(m0);
            if numel(win)==0  % need to create the window
                win=SquareWindow(n0,n0(1)/windowFraction);
                m=single(zeros([n0 nfiles])); % output array
            end;
            %         if pixA==0
            %             error('pixA value is zero.');
            %         end;
            % Window and remove mean. Convert the final images to single to save space.
            m1=m0.*win/(cpe*pixA0^2);  % m1 is in units of electrons
            dose=sum(m1(:))/sum(win(:));
            % remove the mean and store as e / A^2
            m(:,:,findex)=single(m1-win*dose);
            doses(findex)=dose;
        end;
    end;
end;
