classdef parallelGridder < gridderObj
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %   The concrete class that implements the parallel gridder
    %   
    %
    %   Written by Hemant D. Tagare, 11/9/2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   setVol creates a cpuVolGrid and 
        %   sets the FTCAS to the user provided vol
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function volGrid=setVol(obj,vol)
            volGrid=cpuVolFourierGrid();
            volGrid=volGrid.setVolFTCAS(vol);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Extract all images by creating a parallel image 
        %   extractor and using it with volGrid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function imgs=extractImgs(obj,volGrid,angles)
            imgE=parallelImgExtractor();
            imgs=imgE.extractAllImgs(volGrid,angles);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Insert all images by creating a serial image inserter
        %   and using it
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function insertImgs(obj,volGrid,imgs,angles)
%             imgI=parallelImgInserter();
%             imgI.insertAllImgs(volGrid,imgs,angles);
%         end
        function insertImgs(obj,volGrid,imgs,angles,varargin)
            imgI=parallelImgInserter();
            if isempty(varargin)
                ctfs=[];    %No ctfs
            else
                ctfs=[varargin{1,1}];
            end
            imgI.insertAllImgs(volGrid,imgs,angles,ctfs);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Get volume converts the FTCAS back to vol
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function vol=getVol(obj,volGrid)
            vol=volGrid.getVolFTCAS();
        end
    end
end