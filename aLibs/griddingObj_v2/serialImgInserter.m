classdef serialImgInserter 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class makes a concrete image inserter out of the dataResampler
    % abstract class by implmenting the appropriate dataSwapper
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        imgI=[];
    end

    
    methods 
        function insertAllImgs(obj,volGrid,imgs,angles,varargin)
            imgI=imgInserter(volGrid);
            nAngles=size(angles,1);
            ctfs=[varargin{1,1}];
            %Insert all images
            for i=1:nAngles,
                tmpImg=squeeze(imgs(:,:,i));
                if isempty(ctfs)
                    tmpCtf=ones(size(tmpImg));
                else
                    tmpCtf=squeeze(ctfs(:,:,1));
                end
                imgI.insertImg(tmpImg,angles(i,:),tmpCtf);
            end
        end

   end
end