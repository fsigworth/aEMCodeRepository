function x = iswt3(wt,varargin)
%IDWT3 

%   A. Kucukelbir 06-July-2011

% Check arguments.
msg = nargchk(1,2,nargin); %#ok<NCHK>
if ~isempty(msg)
    error('Wavelet:FunctionInput:NbArg',msg)
end
s   = wt.sizeINI;
dec = wt.dec;
Lo  = wt.filters.LoR;
Hi  = wt.filters.HiR;
dwtEXTM = wt.mode;
perFLAG = isequal(dwtEXTM,'per');

lf = zeros(1,3);
for k = 1:3 , lf(k) = length(Lo{k}); end

k = 1;
while k<=length(varargin)
    if ischar(varargin{k})
        word = varargin{k};
        switch word
            case {'a','l','A','L','1'}
                dec{1,1,2}(:) = 0; 
                dec{1,2,1}(:) = 0; dec{1,2,2}(:) = 0;
                dec{2,1,1}(:) = 0; dec{2,1,2}(:) = 0; 
                dec{2,2,1}(:) = 0; dec{2,2,2}(:) = 0;
                
            case {'d','h','D','H','0'}
                dec{1,1,1}(:) = 0; 
                
            otherwise
                if length(word)==3
                    num = ones(1,3);
                    for k = 1:3
                        switch word(k)
                            case {'a','l','A','L','1'}
                            case {'d','h','D','H','0'} , num(k) = 2;
                            otherwise , num(k) = -1; % ERROR
                        end
                    end
                    for n=1:2
                        for j=1:2
                            for k = 1:2
                                if ~isequal([n,j,k],num)
                                    dec{n,j,k}(:) = 0;
                                end
                            end;
                        end;
                    end;
                else
                    error('Wavelet:FunctionInput:ArgVal', ...
                        'Invalid argument value!');
                end
        end
        k = k+1; 
    else
        s = varargin{k};
        k = k+1;
    end
end

% Reconstruction.
perm = [1,3,2];
V = cell(2,2);
for i = 1:2    
    for j = 1:2
        V{j,i} = 0.5* ( swrec1D(dec{i,j,1},Lo{3},perm,perFLAG,s) + ...
                 swrec1D(dec{i,j,2},Hi{3},perm,perFLAG,s) );
    end
end
perm = [2,1,3];
W = cell(1,2);
for i = 1:2
    W{i} = 0.5 * ( swrec1D(V{i,1},Lo{2},perm,perFLAG,s) + ...
        swrec1D(V{i,2},Hi{2},perm,perFLAG,s) );
end

% Last reconstruction.
x = 0.5 * ( swrec1D(W{1},Lo{1},[],perFLAG,s) + swrec1D(W{2},Hi{1},[],perFLAG,s) );

%-----------------------------------------------------------------------%
function X = swrec1D(X,F,perm,perFLAG,s)

if ~isempty(perm)
    X = permute(X,perm);
    s = s(perm);
end
if perFLAG
    lf = length(F);
    lx = size(X,2);
%     lc = lx+lf-1;
%     I = [lx-lf+1:lx , 1:lx , 1:lf];
    nb = lf-1; %fix(lf/2-1);
    idxAdd = 1:nb;
    if nb>lx
        idxAdd = rem(idxAdd,lx);
        idxAdd(idxAdd==0) = lx;
    end
    X = [X X(:,idxAdd,:)];
%     X = X(:,I,:);
end
% sX = size(X);
% if length(sX)<3 , sX(3) = 1; end
% Z = zeros(sX(1),2*sX(2)-1,sX(3));
% Z(:,1:2:end,:) = X;
X = convn(X,F);

% sX = size(X,2);
% F  = floor((sX-s)/2);
% C  = ceil((sX-s)/2);
% X  = X(:,1+F(2):end-C(2),:);

% lenX = size(X,2);
% first = lf; 
% last = lenX-lf+1;
% X = X(:,first:last,:);
% lenX = size(X,2);
% first = 1+floor((lenX-lc)/2);  
% last = first+lc-1;
% X = X(:,first:last,:); 
% 
% last = ceil(lx);
% X = X(:,1:last,:);

X = X(:, lf:lx+lf-1 ,:);

if ~isempty(perm) , X = permute(X,perm); end
%-----------------------------------------------------------------------%
