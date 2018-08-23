classdef cpuMemMapper < memMapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The concrete cpu memory mapper
%
%
%   Copyright Hemant D. Tagare, 2015. Written 11.7.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   This method allocates real or complex zeros depending
        %   on whether varargin{2} is set to 'complex'.
        %   The default is real
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x=allocZeros(obj,varargin)
            
            if nargin==2
                %Allocate real zeros
                x=zeros([varargin{1}]);
            end
            if nargin>2 & strcmp([varargin{2}],'complex')
                %Allocate complex zeros
                x=complex(zeros([varargin{1}]),zeros([varargin{1}]));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Null methods that do not implement anything because
        %   x is in cpu memory already
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x=allocTransfer(obj,x) %Null function
        end
        function x=gatherMem(obj,x)
        end
    end
end