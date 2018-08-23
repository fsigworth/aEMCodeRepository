classdef convKernelObj
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The convolution kernel object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        w=[] %Fred's w-table
        kernelSize=3; %convolution kernel size
        ovr=1024;  %Oversampling rate
        ovrCtr=0; %Center of the table
    end
    
    methods
        function obj=convKernelObj(memMapper)
            %obj.kernelSize=kernelSize
            obj.w=makeKaiserTable(obj.kernelSize,obj.ovr);
            obj.w=obj.w';
            obj.w=memMapper.allocTransfer(obj.w);
           % obj.ovr=ovr;
            obj.ovrCtr=obj.ovr/2+1;
        end
        
        function w1=getConvWeights(obj,i,j,k,iOffset,jOffset,kOffset)
        	iInt=round(i);
            iFrac=floor((i-iInt)*obj.ovr)+obj.ovrCtr;
            jInt=round(j);
            jFrac=floor((j-jInt)*obj.ovr)+obj.ovrCtr;
            kInt=round(k);
            kFrac=floor((k-kInt)*obj.ovr)+obj.ovrCtr;
            nw2=floor(obj.kernelSize/2);
            w1=obj.w(iFrac,iOffset+nw2+1).*obj.w(jFrac,jOffset+nw2+1).*obj.w(kFrac,kOffset+nw2+1);
        end
    end
    
end %End classdef