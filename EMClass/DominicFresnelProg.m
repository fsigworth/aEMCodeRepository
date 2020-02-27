%% 1D step-by-step Fresnel propagation
% Program for step-by-step fresnel propagation from slit in 1d
% clear all;
% close all
width=input('What is the width of the slit? (mm)');
wl=input('What is the wavelength of light (nm)?' )*1e-9;
wdth=1e-3;
origx=[-20:0.01:20]*wdth;
x=origx;
origy=abs(x)<=(width*wdth)/2;
y=origy;
k=2*pi/wl;

n=input('Select number of propagation steps');
dz=input('Select distance of propagation steps (m)');
for S=1:n;
    z1=(S-1)*dz;
    z2=z1+dz;
    [x,y,z1]=fresnel1b(x,y,wl,z1,dz);
    y2=y/max(abs(y));
    scaling1=(wl)*z2/((origx(2)-origx(1))*(length(origx)-1));
    endx=[-((length(origx)-1)/2):((length(origx)-1)/2)]*scaling1;
    h1=((exp(((1i*k)/(2*z2))*origx.^2)));
    yscal=sum(((origy.*abs(h1)).^2)*(origx(2)-origx(1)))/sum((abs(fft(origy.*h1)).^2)*(endx(2)-endx(1)));
    h2=(1/(sqrt(1i*(wl)*z2)))*(exp(((1i*k)/(2*z2))*(endx.^2)));
    fres2=fftshift(fft(origy.*h1))*sqrt(yscal).*h2;
    fres2n=abs(fres2)/max(abs(fres2));
    fresang=angle(fres2);
    
    figure(S)
    subplot(1,3,1)
    plot(x,y2)
    title(['Field distribution after ',num2str(S),' propagation steps of ',num2str(dz),'m']);
    xlabel('x (m)');
    ylabel('Normalised Amplitude');
    subplot(1,3,2)
    plot(x,fres2n)
    title(['Field distribution computed with single step propagation to ',num2str(z2),'m']);
    xlabel('x (m)');
    ylabel('Normalised amplitude');
    subplot(1,3,3)
    error=(fres2n.^2)*(endx(2)-endx(1))-(y2.^2)*(x(2)-x(1));
    plot(endx,error)
    title(['Error between single and ',num2str(S),' step propagation to ',num2str(z2),'m']);
    xlabel('x (m)');
    ylabel('Error');
    sumerror=sum(error);
    percdiff=(abs(((fres2n)./y2)))/max(abs(((fres2n)./y2))); 
    percdiffmean(S)=mean(percdiff);
    CHISQ(S)=sum(((y2-fres2n).^2)./fres2n);
    sumerrorArr(S)=sumerror;
    MEANY(S)=mean(y2);
    MEANF(S)=mean(fres2n);
    hold on
end
figure(n+1)

plot(sumerrorArr);
title('Graph showing accumulation of error with increasing number of propagation steps');
xlabel('Propagation step');
ylabel('Difference in square integral');
figure(n+2)
plot(CHISQ)
hold off
for I=1:25
    d=2*I;
    EVEN(I)=CHISQ(d);
    hold on
end
for I=1:25
    d=2*I-1;
    ODD(I)=CHISQ(d);
    hold on
end
hold off
figure(n+3)
plot(EVEN)
title('X^2 vs even propagation step');
xlabel('Propagation step')
ylabel('X^2')
figure(n+4)
plot(ODD)
title('X^2 vs odd propagation step');
xlabel('Propagation step')
ylabel('X^2')

return




% Function for 1D Fresnel propagation a distance z given spatial field
% distribution, wavelength, propagation distance. For use as operator in
% wave-optics analysis
function [endx,abfres,z2]=fresnel1b(x,y,wl,z1,dz);
if z1~=0
    z2=z1+dz;
else
    z2=dz;
end
k=(2*pi)/(wl);
scaling1=(wl*z2)/((x(2)-x(1))*(length(x)-1));
endx=[-((length(x)-1)/2):((length(x)-1)/2)]*scaling1;
h1=((exp(((1i*k)/(2*z2))*(x.^2)))); h2=(1/(sqrt(1i*(wl)*dz)))*(exp(((1i*k)/(2*dz))*(endx.^2)));
yscal=sum(((y.*abs(h1)).^2)*(x(2)-x(1)))/sum((abs(fft(y.*h1)).^2)*(endx(2)-endx(1)));

fres2=fftshift(fft(y.*h1)).*h2*sqrt(yscal);
abfres=abs(fres2);
end

% Function for 2D Fresnel propagation a distance z given spatial field
% distribution, wavelength, propagation distance. For use as operator in
% wave-optics analysis
function [MM,endx,endy,z2]=fresnel2(M,x,y,wl,z1,dz);
if z1~=0
    z2=z1+dz;
else
    z2=dz;
end
scaling=((wl)*z2)/((x(2)-x(1))*(length(x)-1));
k=(2*pi)/wl;
[X,Y]=meshgrid(x,y);
h2=(exp(((1i*k)/(2*z2))*((X.^2)+(Y .^2))));
endx=[-((length(x)-1)/2):((length(x)-1)/2)]*((wl)*z2/((length(x)-1)*scaling));
endy=[-((length(y)-1)/2):((length(y)-1)/2)]*((wl)*z2/((length(x)-1)*scaling));
[XXX,YYY]=meshgrid(endx,endy);
mult=((exp(1i*k*dz))/(1i*(wl)*dz))*(exp(((1i*k)/(2*dz))*(XXX.^2 + YYY.^2)));
fres2=fftshift(fft2(M.*h2));
yscal=(sum(sum((M.^2)*(x(2)-x(1))*(y(2)-y(1)))))/(sum(sum((fres2.^2)*(endx(2)-endx(1))*(endy(2)-endy(1)))));
fres2r=sqrt(yscal)*fftshift(fft2(M.*h2));
MM=abs(mult.*fres2r); figure(2)
imagesc(endx,endy,MM);
axis image
colormap(hot)
title(['Amplitude distribution in real space at distance ',num2str(z2),'m']); xlabel('x position (m)');ylabel('y position (m)')

end



