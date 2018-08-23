classdef Box < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   This class defines a volume (Box)
    %   The Box maybe contained in a bigger Box called the the container
    %   And the Box may contain another Box
    %   Box is derived from the handle class so the Box is a pointer to the
    %   instance. This is used in the setContainer method to set the
    %   containedBox of the container
    %
    %   Copyright Hemant D. Tagare, 2015. Written 10.27.2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties              %Access to properties is public
        data=[];            %This holds the data         
        boxSize=0;          %Size of the box
        ctr=0;              %Center of the box
        containerBox=[];    %The container of this box
        shiftInContainer=0; %The shift index where this box is contained 
                            %in the container
        indexBeginInContainer=0; %Beginning index of the box in the container
        indexEndInContainer=0;  %End index of the box in the container
        containedBox=[];    %A box contained in this box
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Constructor for Box. If there is an argument, the use the argument
        %   as size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=Box(varargin)
            if nargin>=1
                obj=obj.setSize(varargin{1});
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Utility methods to set and get size and centers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=setSize(obj,boxSize)
            obj.boxSize=boxSize;
            obj.ctr=boxSize/2+1;
        end
        function [boxSize,ctr]=getSize(obj)
            boxSize=obj.boxSize;
            ctr=obj.ctr;
        end
        %Method to set the data
        function obj=setData(obj,x,memMapper)
            obj.data=memMapper.allocTransfer(x); %Assumes that a memory 
                                                 %class is already defined
            obj=obj.setSize(size(x,1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Method to set the container of the box. The method also sets the
        %   contained property of the container to the current Box. The method
        %   relies on the container being a handle object (pointer)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=setContainer(obj,b)
            obj.containerBox=b;
            [containerSize,containerCtr]=b.getSize();
            obj.shiftInContainer=containerCtr-obj.ctr;
            obj.indexBeginInContainer=obj.shiftInContainer+1;
            obj.indexEndInContainer=obj.shiftInContainer+obj.boxSize;
            b.containedBox=obj;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Method duplicates the current Box. It copies over the 
        %   data, sets the sizes and center, but does copy the container
        %   and contained boxes. Used  to duplicate boxes for 
        %   the parallel version of gridding code
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function box1=duplicateBox(obj)
            %Does not copy Container Box
            box1=Box();
            if ~isempty(obj.data)
                box1.data=obj.data;            %This holds the data    
            end
            box1.boxSize=obj.boxSize;          %Size of the box
            box1.ctr=obj.ctr;              %Center of the box
        end
    end
end