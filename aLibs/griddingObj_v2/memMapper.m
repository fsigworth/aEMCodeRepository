classdef memMapper
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Abstract object to create zeros, allocate and transfer
    %   variables, and gather back variables from memory. The memory
    %   maybe cpu or gpu, depending on the concrete class instantiated
    %
    %
    %   Copyright Hemant D. Tagare, 2015. Written 11.7.2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Abstract)
        allocZeros(obj,varargin)
        allocTransfer(obj,x)
        gatherMem(x)
    end
end