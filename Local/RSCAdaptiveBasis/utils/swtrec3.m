function X = swtrec3(wdec,varargin)
%SWTREC3 

%   A. Kucukelbir 07-July-2011

% Check arguments.
nbIn = length(varargin);

% Initialization.
cfs   = wdec.dec;
sizes = wdec.sizes;
level = wdec.level;
cfsFLAG = false;
nbREC_UP = level;

if nbIn>0
    errorFLAG = false;
    type = varargin{1};
    nbChar = length(type);
    if nbChar==2 || nbChar==4
        if isequal(upper(type(1)),'C')
            type(1) = []; nbChar = nbChar-1; cfsFLAG = true;
        else
            errorFLAG = true;
        end
    end
    switch nbChar
        case {1,3}
            num = ones(1,3);
            switch type
                case {'a','l','A','L','1'} , num = Inf;
                case {'d','h','D','H','0'} , num = -num;
                otherwise
                    if nbChar==3
                        for k = 1:3
                            switch type(k)
                                case {'a','l','A','L','1'}
                                case {'d','h','D','H','0'} , num(k) = 2;
                                otherwise , errorFLAG = true;
                            end
                        end
                    else
                        errorFLAG = true;
                    end
            end                        
        otherwise , errorFLAG = true;
    end
    if isequal(num,[1 1 1]) , num = Inf; end
    
    if errorFLAG
        error('Wavelet:FunctionInput:ArgVal', ...
            'Invalid argument value!');
    end
    
    if nbIn>1
        levREC = varargin{2};
        OKval = isnumeric(levREC) && isequal(levREC,fix(levREC)) && ...
            levREC<=level && (levREC>0 || ...
            (levREC==0 && (isequal(num,[1 1 1]) || isinf(num))) );
        if ~OKval
                error('Wavelet:FunctionInput:ArgVal', ...
                    'Invalid argument value!');
        end
    else
        levREC = level;
    end
    if isinf(num)
        First_toDEL = 2+7*(level-levREC);
        for k = First_toDEL:(7*level+1) , cfs{k}(:) = 0; end
        if cfsFLAG
            if isequal(level,levREC) , X = cfs{1}; return; end
            nbREC_UP = level-levREC;
        end
        
    elseif isequal(num,[-1,-1,-1])
        cfs{1}(:) = 0;
        for k = 1:(level-levREC)
            for j = 1:7 , cfs{7*(k-1)+j+1}(:) = 0; end
        end
        if cfsFLAG , X = cfs; end
        
    else
        num(num==2) = 0;
        Idx_toKeep = 1+7*(level-levREC+1) - bin2dec(dec2bin(num)');
        if cfsFLAG
            if (~isequal(num,[1,1,1]) || isequal(level,levREC))
                X = cfs{Idx_toKeep};
                return
            elseif isequal(num,[1,1,1])
                nbREC_UP = level-levREC;
            end
        end
        for k = 1:(7*level+1)
            if k~=Idx_toKeep , cfs{k}(:) = 0; end
        end
        
    end
end


evenoddVal = 0;
idxBeg = 1;
for k=1:nbREC_UP
    idxEnd = idxBeg+7;
    wdec.dec = reshape(cfs(idxBeg:idxEnd),2,2,2);
    X = iswt3(wdec,sizes(k+1,:));
    cfs{idxEnd} = X;
    idxBeg = idxEnd;
    
    % downsample filters.
    for i = 1:size(wdec.filters.LoR,2)
      wdec.filters.LoR{i} = dyaddown(wdec.filters.LoR{i},evenoddVal);
      wdec.filters.HiR{i} = dyaddown(wdec.filters.HiR{i},evenoddVal);
    end
    
end
