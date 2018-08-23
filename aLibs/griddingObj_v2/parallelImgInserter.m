classdef parallelImgInserter 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class makes a concrete image inserter out of the dataResampler
    % abstract class by implmenting the appropriate dataSwapper
    %
    % Parallel image insertion is a bit complicated. The volGrid is
    % duplicated as many times as the number of Matlab workers in
    % volGridCell. Then a subset of images along with each volGridCell 
    % entry is passed to interImgsForParallel which does the serial
    % insertion for each worker. The volGrids are then accumulated in the
    % initial volGrid for the final result.
    %
    % Note that this means that the volGrid is duplicated nWorker times
    % and the images are duplicated once
    %
    %   Written by Hemant D. Tagare, 11/9/2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        imgI=[];
    end

    
    methods 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   This method inserts all the images into the volGrid at
        %   the given angles in parallel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function insertAllImgs(obj,volGrid,imgs,angles,varargin)
            %Decide about ctfs
            if isempty(varargin)
                ctfs=ones(size(imgs));    %No ctfs
            else
                ctfs=[varargin{1,1}];
            end
            %Assume that the parallel pool has started
            % Get the number of workers
            pool=gcp('nocreate');
            nW=pool.NumWorkers;
            %Figure out the image allocation to each worker
            nImgs=size(imgs,3);
            imgPerWorker=floor(nImgs/nW);
            imgBegin=[0:nW-1]*imgPerWorker+1;
            imgEnd=min(imgBegin-1+imgPerWorker,nImgs);
            imgEnd(end)=nImgs;
            %Duplicate the VolGrid
            volGridCell=cell(nW,1);
            for i=1:nW,
                volGridCell{i}=volGrid.duplicate();
            end
            %Insert each allocated image into the volGridCell
            parfor i=1:nW,
                  volGridCell{i}=insertImgsForParallel(volGridCell{i},...
                    imgs(:,:,[imgBegin(i):imgEnd(i)]),...
                    angles([imgBegin(i):imgEnd(i)],:),...
                    ctfs(:,:,[imgBegin(i):imgEnd(i)]));
            end
            %Accumulate the results
            for i=1:nW,
                volGrid.paddedBox.data=volGrid.paddedBox.data+...
                             volGridCell{i}.paddedBox.data;
                volGrid.ctfBox.data=volGrid.ctfBox.data+...   
                            volGridCell{i}.ctfBox.data;
            end
            clear volGridCell;
        end

   end
end