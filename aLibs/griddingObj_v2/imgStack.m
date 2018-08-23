classdef imgStack < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A stack of images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   properties
       imgs=[];
       nImgs=0;
       imgSize=0;
   end
   
   methods
       function obj=imgStack(varargin)
          if nargin==1
              obj=obj.setImgs([varargin{1}]);
          end
       end
       function obj=setImgs(obj,imgs)
           obj.imgs=imgs;
           obj.nImgs=size(imgs,3);
           obj.imgSize=size(imgs,1);
       end
   end
end