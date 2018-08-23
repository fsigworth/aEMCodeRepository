function wt = swt3(x,varargin)
%SWT3

%   A. Kucukelbir 06-July-2011

% Check arguments.
msg = nargchk(2,4,nargin); %#ok<NCHK>
if ~isempty(msg)
    error('Wavelet:FunctionInput:NbArg',msg)
end
LoD = cell(1,3); HiD = cell(1,3); LoR = cell(1,3); HiR = cell(1,3);
argStatus = true;
nextARG = 2;
if ischar(varargin{1})
    [LD,HD,LR,HR] = wfilters(varargin{1}); 
    for k = 1:3
        LoD{k} = LD; HiD{k} = HD; LoR{k} = LR; HiR{k} = HR;
    end
    
elseif isstruct(varargin{1})
    if isfield(varargin{1},'w1') && isfield(varargin{1},'w2') && ...
            isfield(varargin{1},'w3')
        for k = 1:3
            [LoD{k},HiD{k},LoR{k},HiR{k}] = ...
                wfilters(varargin{1}.(['w' int2str(k)]));
        end
    elseif isfield(varargin{1},'LoD') && isfield(varargin{1},'HiD') && ...
           isfield(varargin{1},'LoR') && isfield(varargin{1},'HiR')
        for k = 1:3
            LoD{k} = varargin{1}.LoD{k}; HiD{k} = varargin{1}.HiD{k};
            LoR{k} = varargin{1}.LoR{k}; HiR{k} = varargin{1}.HiR{k};
        end
    else
        argStatus = false;
    end
    
elseif iscell(varargin{1})
    if ischar(varargin{1}{1})
        for k = 1:3
            [LoD{k},HiD{k},LoR{k},HiR{k}] = wfilters(varargin{1}{k});
        end
    elseif iscell(varargin{1})
        Sarg = size(varargin{1});
        if isequal(Sarg,[1 4])
            if ~iscell(varargin{1}{1})
                LoD(1:end) = varargin{1}(1); HiD(1:end) = varargin{1}(2);
                LoR(1:end) = varargin{1}(3); HiR(1:end) = varargin{1}(4);
            else
                LoD = varargin{1}{1}; HiD = varargin{1}{2};
                LoR = varargin{1}{3}; HiR = varargin{1}{4};
            end
        elseif isequal(Sarg,[3 4])
            LoD = varargin{1}(:,1)'; HiD = varargin{1}(:,2)';
            LoR = varargin{1}(:,3)'; HiR = varargin{1}(:,4)';
        else
            argStatus = false;
        end
    end
else
    argStatus = false;
end
if ~argStatus
    error('Wavelet:FunctionInput:ArgVal','Invalid argument value!');
end
sX = size(x);

% Check arguments for Extension.
dwtEXTM = 'per';
for k = nextARG:2:nargin-1
    switch varargin{k}
      case 'mode'  , dwtEXTM = varargin{k+1};
    end
end

% Compute stationary wavelet coefficients.
x = double(x);
dec = cell(2,2,2);
permVect = [];
[a_Lo,d_Hi] = swdec1D(x,LoD{1},HiD{1},permVect,dwtEXTM);

permVect = [2,1,3];
[aa_Lo_Lo,da_Lo_Hi] = swdec1D(a_Lo,LoD{2},HiD{2},permVect,dwtEXTM);
[ad_Hi_Lo,dd_Hi_Hi] = swdec1D(d_Hi,LoD{2},HiD{2},permVect,dwtEXTM);

permVect = [1,3,2];
[dec{1,1,1},dec{1,1,2}] = swdec1D(aa_Lo_Lo,LoD{3},HiD{3},permVect,dwtEXTM);
[dec{2,1,1},dec{2,1,2}] = swdec1D(da_Lo_Hi,LoD{3},HiD{3},permVect,dwtEXTM);
[dec{1,2,1},dec{1,2,2}] = swdec1D(ad_Hi_Lo,LoD{3},HiD{3},permVect,dwtEXTM);
[dec{2,2,1},dec{2,2,2}] = swdec1D(dd_Hi_Hi,LoD{3},HiD{3},permVect,dwtEXTM);

wt.sizeINI = sX;
wt.filters.LoD = LoD;
wt.filters.HiD = HiD;
wt.filters.LoR = LoR;
wt.filters.HiR = HiR;
wt.mode = dwtEXTM;
wt.dec = dec;

%-----------------------------------------------------------------------%
function [L,H] = swdec1D(X,Lo,Hi,perm,dwtEXTM)

if ~isempty(perm) , X = permute(X,perm); end
sX = size(X);
if length(sX)<3 , sX(3) = 1; end

lf = length(Lo);
lx = sX(2);
lc = lx+lf-1;
if lx<lf+1
    nbAdd = lf-lx+1;
    switch dwtEXTM
        case {'sym','symh','symw','asym','asymh','asymw','ppd'}
            Add = zeros(sX(1),nbAdd,sX(3));
            X = [Add , X , Add];
    end
end

switch dwtEXTM
    case 'zpd'             % Zero extension.
        
    case {'sym','symh'}    % Symmetric extension (half-point).
        X = [X(:,lf-1:-1:1,:) , X , X(:,end:-1:end-lf+1,:)];
        
    case 'sp0'             % Smooth extension of order 0.
        X = [X(:,ones(1,lf-1),:) , X , X(:,lx*ones(1,lf-1),:)];
        
    case {'sp1','spd'}     % Smooth extension of order 1.
        Z = zeros(sX(1),sX(2)+ 2*lf-2,sX(3));
        Z(:,lf:lf+lx-1,:) = X;
        last = sX(2)+lf-1;
        for k = 1:lf-1
            Z(:,last+k,:) = 2*Z(:,last+k-1,:)- Z(:,last+k-2,:);
            Z(:,lf-k,:)   = 2*Z(:,lf-k+1,:)- Z(:,lf-k+2,:);
        end
        X = Z; clear Z;
        
    case 'symw'            % Symmetric extension (whole-point).
        X = [X(:,lf:-1:2,:) , X , X(:,end-1:-1:end-lf,:)];
        
    case {'asym','asymh'}  % Antisymmetric extension (half-point).
        X = [-X(:,lf-1:-1:1,:) , X , -X(:,end:-1:end-lf+1,:)];        
        
    case 'asymw'           % Antisymmetric extension (whole-point).
        X = [-X(:,lf:-1:2,:) , X , -X(:,end-1:-1:end-lf,:)];

    case 'rndu'            % Uniformly randomized extension.
        X = [randn(sX(1),lf-1,sX(3)) , X , randn(sX(1),lf-1,sX(3))];        
                        
    case 'rndn'            % Normally randomized extension.
        X = [randn(sX(1),lf-1,sX(3)) , X , randn(sX(1),lf-1,sX(3))];        
                
    case 'ppd'             % Periodized extension (1).
        X = [X(:,end-lf+2:end,:) , X , X(:,1:lf-1,:)];
        
    case 'per'             % Periodized extension (2).
        if rem(lx,2) , X = [X , X(:,end,:)]; lx = lx + 1; end
        I = [lx-lf+1:lx , 1:lx , 1:lf];
        if lx<lf
            I = mod(I,lx);
            I(I==0) = lx;
        end
        X = X(:,I,:);
end

L = convn(X,Lo);
H = convn(X,Hi);
clear X
switch dwtEXTM
    case 'zpd'
    otherwise
        lenL = size(L,2);
        first = lf; last = lenL-lf+1;
        L = L(:,first:last,:); H = H(:,first:last,:);
        lenL = size(L,2);
        first = 1+floor((lenL-lc)/2);  last = first+lc-1;
        L = L(:,first:last,:); H = H(:,first:last,:);
end
% L = L(:,2:2:end,:);
% H = H(:,2:2:end,:);
if isequal(dwtEXTM,'per')
    last = ceil(lx);
    L = L(:,1:last,:);
    H = H(:,1:last,:);
end

if ~isempty(perm)
    L = permute(L,perm);
    H = permute(H,perm);
end
%-----------------------------------------------------------------------%
