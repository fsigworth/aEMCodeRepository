function ctfInt=ctfInterp(ctf,volBox,interpBox)
            tmp=complex(ctf,zeros(size(ctf)));
            tmp=fftshift(real(ifft2(fftshift(tmp))));
            imgInterp=zeros([1 1]*interpBox.boxSize);
            iBegin=volBox.indexBeginInContainer;
            iEnd=volBox.indexEndInContainer;
            imgInterp(iBegin:iEnd,iBegin:iEnd)=tmp;
            ctfInt=real(fftshift(fft2(fftshift(imgInterp))));
end