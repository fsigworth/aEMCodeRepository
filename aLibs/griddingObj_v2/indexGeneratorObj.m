classdef indexGeneratorObj
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This object does the translation from the plane coordinates to the 
%   volume coordinates. The object is initialized with VolGrid. The
%   offsetAddrs method takes plane i and j coordinates of grid points
%   incremented at stride and starting from the gridStart coordinates
%   and returns the plane and volume coordinates. The stride and offset are
%   used to vectorize the dataResampling
%
%   The object uses Fred Sigworth's EulerMatrix function.
%
%   Written by Hemant D. Tagare, 11/9/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   properties
        volGrid=[];
        %Properties to hold variables for index
        % calculation
        %imgIndex=[];
        imgCtr=0;
        volIndexLims=zeros(3,2);
        volCtr=0;
        %Properties to hold variables for the 
        %orientation of the plane
        E=[];       %Euler matrix
        %Convolution Kernel
        kernelSize=0;
        ovr=0;
        ovrCtr=0;
   end %End properties
    
   methods
       
       function obj=indexGeneratorObj(volGrid)
            obj.volGrid=volGrid;
            obj.kernelSize=volGrid.convKernel.kernelSize;
            obj.ovr=volGrid.convKernel.ovr;
            obj.ovrCtr=volGrid.convKernel.ovrCtr;
            %Allocate memory for interpolated image
            interpSize=volGrid.interpBox.boxSize;
            %Set the limits for index mapping
            obj.volIndexLims(:,1)=volGrid.interpBox.indexBeginInContainer;
            obj.volIndexLims(:,2)=volGrid.interpBox.indexEndInContainer;
            %Set the centers
            obj.imgCtr=volGrid.interpBox.ctr;
            obj.volCtr=volGrid.paddedBox.ctr;
       end
       
       function obj=setAngles(obj,angles)
            %Set the Euler matrix
            obj.E=inv(EulerMatrix(angles));
       end
       
       function [planeAddr,volAddr]=offsetAddrs(obj,stride,gridStart)
                %Plane addrs
                [iPlane,tmp]=ndgrid(gridStart(1):stride:obj.volGrid.interpBox.boxSize);
                iPlane=reshape(iPlane,[numel(iPlane) 1]);
                planeAddr=iPlane;
                iPlane=iPlane-obj.volGrid.interpBox.ctr;
                [tmp,jPlane]=ndgrid(gridStart(2):stride:obj.volGrid.interpBox.boxSize);
                jPlane=reshape(jPlane,[numel(jPlane) 1]);
                planeAddr=[planeAddr jPlane];
                jPlane=jPlane-obj.volGrid.interpBox.ctr;
                
                %Volume addrs
                iVol=obj.E(1,1)*iPlane+obj.E(1,2)*jPlane+obj.volGrid.paddedBox.ctr;
                iVol=min(max(iVol,obj.volIndexLims(1,1)),obj.volIndexLims(1,2));  % prevent coords from going out of bounds
                
                jVol=obj.E(2,1)*iPlane+obj.E(2,2)*jPlane+obj.volGrid.paddedBox.ctr;
                jVol=min(max(jVol,obj.volIndexLims(1,1)),obj.volIndexLims(1,2));  % prevent coords from going out of bounds
                
                kVol=obj.E(3,1)*iPlane+obj.E(3,2)*jPlane+obj.volGrid.paddedBox.ctr;
                kVol=min(max(kVol,obj.volIndexLims(1,1)),obj.volIndexLims(1,2));  % prevent coords from going out of bounds
                
                volAddr=[iVol jVol kVol];
        end
        
   end%End methods
end%End classdef