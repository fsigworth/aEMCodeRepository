classdef imgExtractor < dataResampler
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class makes a concrete image extractor out of the dataResampler
    % abstract class by implmenting the appropriate dataSwapper
    %
    %   Written by Hemant D. Tagare, 11/9/2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        function img=dataSwapper(obj,volIndexOffset,imgIndex,img,wt,ctf) 
            img(imgIndex)=img(imgIndex)+wt.*(obj.volGrid.paddedBox.data(volIndexOffset)); %Ignore ctf
        end
        
        function obj=imgExtractor(volGrid)
            obj=obj.setGrid(volGrid);
        end
        
        function img=extractImg(obj,img,angles)
            obj=obj.setAngles(angles);
            %Extract the image
            img=obj.volGrid.memMapper.allocZeros([1 1]*obj.volGrid.interpBox.boxSize);
            img=obj.resampleVolImg(img,angles,[]); %No ctfs
            %Resample it
            %Reshape the image and convert back to spatial domain
            img=reshape(img,[obj.volGrid.interpBox.boxSize obj.volGrid.interpBox.boxSize]);
            img=convertFromCAS(img);
            img=fftshift(real(ifft2(fftshift(img))));
            iBegin=obj.volGrid.volBox.indexBeginInContainer;
            iEnd=obj.volGrid.volBox.indexEndInContainer;
            img=img(iBegin:iEnd,iBegin:iEnd);
        end

   end
end