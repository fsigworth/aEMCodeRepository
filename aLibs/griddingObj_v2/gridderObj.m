classdef gridderObj
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   The abstract gridder class. It defines the interface
    %   for concrete gridders
    %
    %
    %
    %   Written by Hemant D. Tagare, 11/9/2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Abstract)
        setVol;         % Set the volume grid
        insertImgs;     % Insert all images
        extractImgs;    % Extract all images
        getVol;         % Get the vol from volGrid
    end
end