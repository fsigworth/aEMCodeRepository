function avgShift=bayesShiftMatch(img1,img2,maxShift,sigma)
    perf=[];
    shiftArray=[];
    for shift=1:maxShift,
        img1=circshift(img1,1,1);
        perf=[perf sum((img1(:)-img2(:)).^2)];
        shiftArray=[shiftArray shift];
    end
    perf=perf-max(perf);
    perf=exp(-perf/(2*sigma^2));
    perf=perf/sum(perf);
    avgShift=sum(shiftArray.*perf);
    


