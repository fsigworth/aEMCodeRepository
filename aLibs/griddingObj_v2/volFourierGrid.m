classdef volFourierGrid < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   The abstract base class to hold the volume Fourier grid
    %   The properties and some methods of the class are 
    %   abstract and have to implemented for specific 
    %   execution - serial, parfor, and gpu
    %
    %   This object holds the CAS Fourier Transform of a volume 
    %   This is an objectified version of Fred Sigworth's gridding code
    %
    %
    %   Copyright Hemant D. Tagare, 2015. Written 11.7.2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Abstract)
        memMapper;              %   This holds the memory mapper object. The memory mapper
                                %   allocates memory and transfers data back
                                %   and forth. Where the memory is allocated
                                %   depends on whether the class derived from
                                %   this class is for a cpu or a  gpu version 
    end
    
    properties      
        interpolationFactor=2;  %   Hardcoded: Interpolate Fourier Transform by 2
        gridPadSize=8;          %   Hardcoded: Eight voxels on each side
        origVolSize=0;
        
        convKernel=[];          %   Convolution kernel Object
        
        volBox=[];              %   Box of original volume size
        interpBox=[];           %   Box of interpolated volume size
        paddedBox=[];           %   Box of interpolated + padded size
        ctfBox=[];              %   Box to hold the CTF weights
    end
    
    methods
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % This method upsamples the spatial domain volume and 
         % converts it into the CAS form in the fourier domain
         % vol is required to be 3d and have equal sides with even
         % number of voxels per side.
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
         function obj=setVolFTCAS(obj,vol)

             %  Test whether the volume is a cube with equal sides that are
             %  even
                if numel(size(vol)) ~= 3
                    disp('****** Error (volFourierGrid): volume is not 3d');
                    return;
                end
                if prod(size(vol)) ~= mean(size(vol))^3
                     disp('****** Error (volFourierGrid): volume does not have equal sides');
                    return;
                end
                if mod(size(vol,1),2) ~= 0
                    disp('****** Error (volFourierGrid): volume is not even');
                    return;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Create Boxes and set up the containment relations
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %   Volume box
                obj.volBox=Box();
                obj.volBox=obj.volBox.setSize(size(vol,1));
                
                %   Create the interpolated Box and set is as a container for
                %   volBox
                obj.interpBox=Box();
                obj.interpBox.setSize(obj.interpolationFactor*obj.volBox.boxSize); % Volume size after interpolation
                obj.volBox.setContainer(obj.interpBox);     %Set the container of the volBox to the interpolation Box;
                
                %   Create the padded Box and set is as a container for
                %   interpBox
                obj.paddedBox=Box();
                obj.paddedBox=obj.paddedBox.setSize(obj.interpBox.boxSize+2*obj.gridPadSize); % Volume size after interpolation
                obj.interpBox.setContainer(obj.paddedBox);     %Set the container of the interpBox to paddedBox;
                
                %   Create the CTF wt Box
                obj.ctfBox=Box();
                obj.ctfBox=obj.ctfBox.setSize(obj.interpBox.boxSize+2*obj.gridPadSize); % Volume size after interpolation
                obj.ctfBox=obj.ctfBox.setData(zeros([1 1 1]*obj.paddedBox.boxSize),obj.memMapper);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Process the volume by 
                % 1) Inserting it into the intepolatedBox
                % 2) Precompensating
                % 3) Interpolating the Fourier Transform
                % 4) Inserting it into the padded box
                % 5) Fuzzymasking the inserted Fourier Transform
                % 6) Converting to CAS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % 1) Inserting into the intepolatedBox
                interpSize=obj.interpBox.boxSize;
                volInterp=zeros(interpSize,interpSize,interpSize);
                indexB=obj.volBox.indexBeginInContainer;
                indexE=obj.volBox.indexEndInContainer;
                volInterp(indexB:indexE,indexB:indexE,indexB:indexE)=vol;
                % 2) Precompensate
                obj.convKernel=convKernelObj(obj.memMapper); %Conv kernel object has the kernel size
                preComp=getPreComp(interpSize, obj.convKernel.kernelSize);
                volInterp=volInterp.*reshape(kron(preComp,kron(preComp,preComp)),...
                    interpSize,interpSize,interpSize);
                % 3) Interpolated Fourier Transform
                interpFT=fftshift(fftn(fftshift(volInterp))); %Interpolate Fourier Transform
                % 4) Insert in a volume of the size of paddedBox
                paddedSize=obj.paddedBox.boxSize;
                paddedFT=complex(zeros(paddedSize,paddedSize,paddedSize),zeros(paddedSize,paddedSize,paddedSize));
                indexB=obj.interpBox.indexBeginInContainer;
                indexE=obj.interpBox.indexEndInContainer;
                paddedFT(indexB:indexE,indexB:indexE,indexB:indexE)=interpFT; %insert the interpolated transform
                % 5) Fuzzy mask the outside 
                paddedFT=paddedFT.*fuzzymask(paddedSize,3,interpSize/2+2,3,...
                                [1 1 1]*obj.paddedBox.ctr); % mask outside frequencies.
                % 6) Convert to CAS and store in paddedBox
                obj.paddedBox=obj.paddedBox.setData(convertToCAS(paddedFT),obj.memMapper);

         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %  Weiner filter the paddedBox data using SNR and the ctfBox data
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function obj=weinerFilt(obj,snr)
             if snr < 1e-12         %   Essentially snr is zero
                 snr= 1e-12;        %   Set it to a small value to avoid
                                    %   division overflow
             end
             obj.paddedBox.data=obj.paddedBox.data./(obj.ctfBox.data+snr);
         end

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %  This method converts the upsampled + padded CAS back into 
         %  the downsampled spatial domain volume. It is the transpose
         %  (adjoint) of the setVolFTCAS operator. The operations are 
         %  executed in reverse order of setVolFTCAS  
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function vol=getVolFTCAS(obj)        
                  % 6) Convert from CAS
                  vol=convertFromCAS(obj.paddedBox.data);
                  % 5) Apply fuzzy mask
                  paddedSize=obj.paddedBox.boxSize;
                  interpSize=obj.interpBox.boxSize;
                  vol=vol.*fuzzymask(paddedSize,3,interpSize/2+2,3,...
                                [1 1 1]*obj.paddedBox.ctr); % mask outside frequencies.
                  % 4) Extract the interpolated part
                  iBegin=obj.interpBox.indexBeginInContainer;
                  iEnd=obj.interpBox.indexEndInContainer;
                  vol=vol(iBegin:iEnd,iBegin:iEnd,iBegin:iEnd);
                  % 3) inverse FFT the interpolated transform
                  vol=fftshift(real(ifftn(fftshift(vol))));
                  % 2) Precompensate
                  convKernel=convKernelObj(obj.memMapper); %Conv kernel object has the kernel size
                  preComp=getPreComp(interpSize, convKernel.kernelSize);
                  vol=vol.*reshape(kron(preComp,kron(preComp,preComp)),...
                    interpSize,interpSize,interpSize);
                  % 1) Extract out the vol at original size
                  iBegin=obj.volBox.indexBeginInContainer;
                  iEnd=obj.volBox.indexEndInContainer;
                  vol=vol(iBegin:iEnd,iBegin:iEnd,iBegin:iEnd);
                  vol=obj.memMapper.gatherMem(vol);
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %  This method duplicates volGrid into the current grid
         %  It calls duplicateBox to duplicate the boxes and then 
         %  sets the containers. It is used to duplicate volGrids in 
         %  parallel code.
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function obj=duplicateGrid(obj,volGrid)
                obj.interpolationFactor=volGrid.interpolationFactor;  
                obj.gridPadSize=volGrid.gridPadSize;       
                obj.origVolSize=volGrid.origVolSize;

                obj.convKernel=volGrid.convKernel;  

                obj.volBox=volGrid.volBox.duplicateBox();           
                obj.interpBox=volGrid.interpBox.duplicateBox();           
                obj.paddedBox=volGrid.paddedBox.duplicateBox();
                obj.ctfBox=volGrid.ctfBox.duplicateBox();
                %Set containers
                obj.volBox.setContainer(obj.interpBox);
                obj.interpBox.setContainer(obj.paddedBox);   
         end
             
    end
end

