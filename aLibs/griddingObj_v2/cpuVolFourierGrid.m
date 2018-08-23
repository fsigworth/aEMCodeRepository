classdef cpuVolFourierGrid < volFourierGrid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A concrete implementation of the abstract volFourierGrid class
    % This holds everything in main memory. The memory mapper is 
    % set to the cpuMemMapper
    %
    %
    %
    %   Copyright Hemant D. Tagare, 2015. Written 11.7.2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties
        memMapper=[];
    end
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Set memory mapper to the cpu Memory mapper
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=cpuVolFourierGrid(obj)
            obj.memMapper=cpuMemMapper();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Duplication for parallelization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj1=duplicate(obj)
            obj1=cpuVolFourierGrid();
            obj1=obj1.duplicateGrid(obj);
        end
    end
end

