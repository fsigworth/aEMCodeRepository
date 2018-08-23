classdef imgInserter < dataResampler
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class makes a concrete image inserter out of the dataResampler
    % abstract class by implmenting the appropriate dataSwapper
    %
    %
    %   Written by Hemant D. Tagare, 11/9/2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        function img=dataSwapper(obj,volIndexOffset,imgIndex,img,wt,ctf)
            obj.volGrid.paddedBox.data(volIndexOffset) = ...
                            obj.volGrid.paddedBox.data(volIndexOffset) + wt.*img(imgIndex).*ctf(imgIndex); %Add the data
            obj.volGrid.ctfBox.data(volIndexOffset) =obj.volGrid.ctfBox.data(volIndexOffset) +wt.*ctf(imgIndex).*ctf(imgIndex);
        end
        
        function obj=imgInserter(volGrid)
            obj=obj.setGrid(volGrid);
        end
        
        function insertImg(obj,img,angles,ctf)
            %Interpolate the image
            imgInterp=imgInterpToCAS(img,obj.volGrid.volBox,obj.volGrid.interpBox);
            %Interpolate the ctf
            ctfInt=ctfInterp(ctf,obj.volGrid.volBox,obj.volGrid.interpBox);
            %Insert it
            obj.resampleVolImg(imgInterp,angles,ctfInt);
        end

   end
end