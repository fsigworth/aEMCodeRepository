function [shift,probab]=bayes_shift_est(model,data,sigma,prRange)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Estimates the 1d circ shift of the model in data %
%   Assume a uniform prior of shift -10 to + 10      %
%   model is an array containing the model           %
%   data is an array containing  the data            %
%   Model and data array have to be same size        %
%   Sigma is the noise stdev                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Prior
    prProbab=zeros(size(model));
    prProbab(1:prRange)=1;
    prProbab(end-prRange+1:end)=1;
    prProbab=prProbab/sum(prProbab);
    %Calculate norms
    modelNorm=norm(model,2)^2;
    dataNorm=norm(data,2)^2;
    corr=real(ifft(fft(data).*conj(fft(model))));
    %Calculate probability
    probab=exp(-(modelNorm-2*corr+dataNorm)/(2*sigma^2)).*prProbab;
    probab=probab/sum(probab);
    %Average over a circle
    theta=2*pi*((1:numel(probab))'-1)/numel(probab);
    sinTheta=sin(theta);
    cosTheta=cos(theta);
    avgY=sum(probab.*sinTheta);
    avgX=sum(probab.*cosTheta);
    meanTheta=atan2(avgY,avgX);
    shift=numel(probab)*meanTheta/(2*pi)+1;
