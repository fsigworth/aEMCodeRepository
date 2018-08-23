function maxIndex=shiftMatch(img1,img2,maxShift)
    perf=[];
    for shift=1:maxShift,
        img1=circshift(img1,1,1);
        perf=[perf sum(img1(:).*img2(:))];
    end
    [maxVal,maxIndex]=max(perf);


