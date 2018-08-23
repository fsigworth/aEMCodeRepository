function maxIndex=frcMatch(img1,img2,maxShift)
    perf=[];
    imgSize=size(img1,1);
    rings=makeRings(imgSize,12);
    t=0.01;
coeffs=zeros(12,maxShift);
for shift=1:maxShift,
        img1=circshift(img1,1,1);
        coeff=frc(img1,img2,rings);
        coeffs(:,shift)=coeff;
        tmp= (coeff>=t);
        perf=[perf sum((coeff-t).*tmp)];
end
    subplot(2,1,1);
    plot(coeffs);
    legend(num2str((1:maxShift)'));
    subplot(2,1,2);
    plot(perf);
pause
%     figure(3);
%     hold on;
%     plot(perf);
%     clear max;
[maxVal,maxIndex]=max(perf);


