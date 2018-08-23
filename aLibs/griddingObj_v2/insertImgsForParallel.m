function volGrid=insertImgsForParallel(volGrid,imgs,angles,ctfs)
            imgI=imgInserter(volGrid);
            nAngles=size(angles,1);
            %Insert all images
            for i=1:nAngles,
                imgI.insertImg(squeeze(imgs(:,:,i)),angles(i,:),squeeze(ctfs(:,:,i)));
            end
end