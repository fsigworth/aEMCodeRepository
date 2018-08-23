function [newModel,probabs]=EM_shift_est(model,data,sigma,prRange)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Estimates the 1d circ shift of the model in data %
%   Assume a uniform prior of shift -10 to + 10      %
%   model is an array containing the model           %
%   data is an array containing  the data            %
%   Model and data array have to be same size        %
%   Sigma is the noise stdev                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nx=100;
% sh=10;
% prRange=20;
% sigma=1;
% xs=(-nx/2:nx/2-1)';
% sigmaM=2;
% m0=exp(-xs.^2/(2*sigmaM^2));
% model=newModel;
% 
% 
% model=m0;
% data=circshift(m0,[sh 0])+sigma*randn(nx,1);
% 
% 
    prProbab=zeros(size(model));
    prProbab(1:prRange)=1;
    prProbab(end-prRange+1:end)=1;
    prProbab=prProbab/sum(prProbab);

probabs=zeros(size(data));
nSet=size(data,2);
newModel=zeros(size(model));
for iSet=1:nSet
    %Prior
    %Calculate norms
    modelNorm=norm(model,2)^2;
    dataNorm=norm(data(:,iSet),2)^2;

    corr=real(ifft(fft(data(:,iSet)).*conj(fft(model))));
    %Calculate probability
    probab=exp(-(modelNorm-2*corr+dataNorm)/(2*sigma^2)).*prProbab;
    probab=probab/sum(probab);
    probabs(:,iSet)=probab;
    newModel=newModel+real(ifft(fft(data(:,iSet)).*conj(fft(probab))));
end;
newModel=newModel/nSet;

%     subplot(2,1,1);
%     plot(probab);
%     subplot(2,1,2);
%     plot([m0 newModel data])
%     
%     return
    
%     
%     %Average over a circle
%     theta=2*pi*((1:numel(probab))'-1)/numel(probab);
%     sinTheta=sin(theta);
%     cosTheta=cos(theta);
%     avgY=sum(probab.*sinTheta);
%     avgX=sum(probab.*cosTheta);
%     meanTheta=atan2(avgY,avgX);
%     shift=numel(probab)*meanTheta/(2*pi)+1;
