function obj=gridder(type,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Essentially a factory that returns the serial or the 
%   parallel gridder. An optional argument for the 
%   parallel gridder sets the number of MATLAB workers
%   This function also sets the path to EMBASE
%
%   Written by Hemant D. Tagare, 11/9/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %path(path,'..\EMBase');
    if strcmp(type,'serial')
        obj=serialGridder();
    end
    if strcmp(type,'parallel')
        %Start the parapool
        pool=gcp('nocreate');
        if isempty(pool) & ~isempty(varargin)
            nW=[varargin{1}];
            parpool(nW);
        end
        if isempty(pool) & isempty(varargin)
            parpool;
        end
        %Return the parallel gridder
        obj=parallelGridder();
    end
            
end