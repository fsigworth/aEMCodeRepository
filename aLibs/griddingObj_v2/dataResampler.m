classdef dataResampler
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This abstract class implements the data resampling 
    % This class captures the idea that inserting or extracting
    % an image into a volume requires the same set of indices. The main
    % difference between the two is the direction in which the data is
    % swapped: From image to volume or backwards. The direction can be 
    % encoded in the Abstract method Data Swapper. 
    %
    % This object is to be used by initializing with a volGrid, followed
    % by setting the plane angle, followed by resampleVolImg. The last
    % uses the initialized volGrid and img as the volume and the plane.
    % 
    %
    %   Written by Hemant D. Tagare, 11/9/2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        volGrid=[];
        imgs=[];
        angles=[];

        %Index generator Obj
        indexGen=[];
    end
    
    methods (Abstract)
        dataSwapper;
    end
    
    methods (Access=private)
       function img=resampleVolImgAtOffset(obj,img,stride,gridStart,ctf)
            %Get the addresses
            [imgAddr,volAddr]=obj.indexGen.offsetAddrs(stride,gridStart);
            %Get sizes 
            interpSize=obj.volGrid.interpBox.boxSize;
            paddedSize=obj.volGrid.paddedBox.boxSize;
           %Generate 1d indices
            imgIndex=imgAddr(:,1)+interpSize*(imgAddr(:,2)-1);
            volIndex=round(volAddr(:,1))+paddedSize*((round(volAddr(:,2))-1)+paddedSize*(round(volAddr(:,3))-1)); 
            %Extract the image
            for iOffset=-1:1,
                for jOffset=-1:1,
                    for kOffset=-1:1,
                        wt=obj.volGrid.convKernel.getConvWeights(volAddr(:,1),volAddr(:,2),volAddr(:,3),iOffset,jOffset,kOffset);
                        volIndexOffset=volIndex+iOffset+paddedSize*((jOffset)+paddedSize*(kOffset));
                        img=obj.dataSwapper(volIndexOffset,imgIndex,img,wt,ctf);
                       
                    end
                end
            end
       end
       
  
   end
    
    methods 
        
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % At initialization volGrid may be used as an 
           % argument
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          function obj=dataResampler(varargin)
                if nargin >=1
                    obj=obj.setGrid(varargin);
                end
          end
          
           function obj=setGrid(obj,volGrid)
            %Store the volGrid,imgs and object
            obj.volGrid=volGrid;
            %Index generator
            obj.indexGen=indexGeneratorObj(volGrid);
          end
        
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Set the plane insertion angle
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          function obj=setAngles(obj,angles)
             obj.angles=angles;
             obj.indexGen=obj.indexGen.setAngles(angles);
          end

   
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %  Resample the volume or the image depending on 
         % the dataswapper
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function img=resampleVolImg(obj,img,angles,ctf)
            %Assume that the volume grid has been set
            %Return img, volGrid is updated automatically
            obj=obj.setAngles(angles);
            %Image is extracted in several part
            nParts=obj.volGrid.convKernel.kernelSize+1;
            for iStart=1:nParts,
                for jStart=1:nParts,
                    img=obj.resampleVolImgAtOffset(img,nParts,[iStart jStart],ctf);  
                end
            end
        end
   end
end