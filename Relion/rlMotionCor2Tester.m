% rlMotionCor2Tester
% Create fake movies for characterization of MotionCor2 processing

signalName='/Users/fred/EMWork/Hideki/170814/KvLipo125_3w10/Merged/sq05_1_0001_Aug14_14.58.28mz.tif';

n=[1 1]*7680;
s0=ReadEMFile(signalName);
s1=GaussFilt(s0,.066)+0.03*s0; % 20 A filter
s1=Downsample(s1,n);
signalLim=5*std(s1(:));
s1=max(-signalLim,min(signalLim,s1));
%%
% noiseSigma=2;
signalAmp=.25;
frameMean=.25;  % had better be larger than signalLim, about 0.2
nFrames=10;
frameSize=size(s1);

driftSigma=0;
driftVelocity=[2 2];

shifts=zeros(nFrames,2);
sOut=zeros([frameSize nFrames],'uint8');
fsignal=fftn(max(0,signalAmp*s1));
for i=1:nFrames
    if any(shifts(i,:))
        p=FourierShift(frameSize,shifts(i,:));
    else
        p=1;
    end;
    s2=real(ifftn(fsignal.*p));
    sOut(:,:,i)=poissrnd(s2+frameMean,frameSize);
    if i<nFrames  % create a random-walk drift for next time
        shifts(i+1,:)=shifts(i,:)+randn(1,2)*driftSigma+driftVelocity;
    end;
end;

disp('Getting a name for the output movie');
[outName,outPath]=uiputfile('*.mrc');
if isnumeric(outName)
    return
end;

disp(['Writing ' AddSlash(outPath) outName]);
cd(outPath);
WriteMRC(sOut,1,outName,0);
disp('display');
imags(GaussFilt(sum(single(sOut),3),.1));
