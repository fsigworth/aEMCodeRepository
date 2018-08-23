classdef serialImgExtractor 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class extracts all images via a serial for loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        imgE=[] %image extractor
    end
    
    methods 
        
        function imgs=extractAllImgs(obj,volGrid,angles)
            imgE=imgExtractor(volGrid);
            nAngles=size(angles,1);
            %Allocate the memory for images on the device
            imgs=volGrid.memMapper.allocZeros([[1 1]*volGrid.volBox.boxSize nAngles]);
            %Extract each image
            for i=1:nAngles,
                imgs(:,:,i)=imgE.extractImg(squeeze(imgs(:,:,i)),angles(i,:));
            end
            %Gather the images
            imgs=volGrid.memMapper.gatherMem(imgs);
        end

   end
end