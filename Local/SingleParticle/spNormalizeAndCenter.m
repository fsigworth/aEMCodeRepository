function [centeredStack,sp,shifts]=spNormalizeAndCenter(...
    stack,pixA,maskRadius,nRounds,displayOn)
% function [centeredStack,sp,shifts]=spNormalizeAndCenter(...
%                           stack,pixA[,maskRadius,nRounds,displayOn])
% Given a stack of particle images, call Normalize them and center them
% against their mean using nRounds iterations (default=3).
% Also create the stack parameters structure using spCreateStackParams.

filtSD=2;  % std (in pixels) of the Gaussian that filters the trans alignment

if nargin<3
    maskRadius=0.48;
end;
if nargin<4
    nRounds=3;
end;

if displayOn
    clf;
    SetGrayscale;
end;


[n, ny, nim]=size(stack);
stack=NormalizeImages(stack,1,0);  % Normalize but don't mask.  Use these for centering.

transSD=n/4;
fc=.133/filtSD;  % Gauss filter frequency
ccprior=-Radius(n).^2/(2*transSD^2);
rmask=fuzzymask(n,2,maskRadius*n,maskRadius*0.1*n);  % basic mask
for i=1:nRounds  % 3 rounds of centering
    ref=CircSymmetrize(mean(stack,3)).*rmask;
    ref=GaussFilt(ref,fc);
    ref=ref/sqrt(ref(:)'*ref(:));  % normalize
    
    if displayOn
        subplot(221);
        imacs(ref);
        title('Normalize and Center');
        subplot(2,2,2);
        imacs(ref);
        title(['Iteration ' num2str(i)]);
        drawnow;
    end;
    [stack, shifts, amps]=TransAlignInt(stack,ref,ccprior);
    %     sum(ok)
    if displayOn
        subplot(2,2,3);
        plot(shifts(1,:),shifts(2,:),'k.');
        axis([-n/2 n/2 -n/2 n/2]);
        title('Shifts');
        subplot(2,2,4);
        hist(amps,100);
        xlabel('Amplitudes');
        drawnow;
    end;
end;
centeredStack=stack;
for i=1:nim
    centeredStack(:,:,i)=stack(:,:,i).*rmask;
end;
sp=spCreateStackParams(centeredStack,pixA);
sp.amp=single(amps);


