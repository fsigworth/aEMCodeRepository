function wdec = swtdec3(X,level,varargin)
%swtdec3 

%   A. Kucukelbir 07-July-2011

% Check arguments.
nbIn = nargin;
LoD = cell(1,3); HiD = cell(1,3); LoR = cell(1,3); HiR = cell(1,3);
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
        error('Wavelet:FunctionInput:ArgVal','Invalid argument value!');
    end
    
elseif iscell(varargin{1})
    if ischar(varargin{1}{1})
        for k = 1:3
            [LoD{k},HiD{k},LoR{k},HiR{k}] = wfilters(varargin{1}{k});
        end
    else
        LoD(1:end) = varargin{1}(1); HiD(1:end) = varargin{1}(2);
        LoR(1:end) = varargin{1}(3); HiR(1:end) = varargin{1}(4);
    end
else
    error('Wavelet:FunctionInput:ArgVal','Invalid argument value!');
end

% Check arguments for Extension.
dwtEXTM = 'per';
for k = 2:2:nbIn-2
    switch varargin{k}
      case 'mode'  , dwtEXTM = varargin{k+1};
    end
end

% Initialization.
evenoddVal = 1;
if isempty(X) , wdec = {}; return; end
sizes = zeros(level+1,3);
sizes(level+1,1:3) = size(X);
for k=1:level
    wdec = swt3(X,{LoD,HiD,LoR,HiR},'mode',dwtEXTM);
    X = wdec.dec{1,1,1};
    if length(size(X))>2
        sizes(level+1-k,1:3) = size(X);
    else
        sizes(level+1-k,1:3) = ceil(sizes(level+2-k,1:3)/2);
    end
    wdec.dec = reshape(wdec.dec,8,1,1);
    if k>1
        cfs(1) = [];
        cfs = cat(1,wdec.dec,cfs);
    else
        cfs = wdec.dec;
    end
    
    % upsample filters.
    for i = 1:size(LoD,2)
      LoD{i} = dyadup(LoD{i},evenoddVal);
      HiD{i} = dyadup(HiD{i},evenoddVal);
      LoR{i} = dyadup(LoR{i},evenoddVal);
      HiR{i} = dyadup(HiR{i},evenoddVal);
    end
    
end
wdec.sizeINI = sizes(end,:);
wdec.level = level;
wdec.dec   = cfs;
wdec.sizes = sizes;
wdec = orderfields(wdec,{'sizeINI','level','filters','mode','dec','sizes'});

