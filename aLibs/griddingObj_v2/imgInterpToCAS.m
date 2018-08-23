function imgInterp=imgInterpToCAS(img,volBox,interpBox)
 
            imgInterp=zeros([1 1]*interpBox.boxSize);
            iBegin=volBox.indexBeginInContainer;
            iEnd=volBox.indexEndInContainer;
            imgInterp(iBegin:iEnd,iBegin:iEnd)=img;
            imgInterp=fftshift(fft2(fftshift(imgInterp)));
            imgInterp=convertToCAS(imgInterp);
end