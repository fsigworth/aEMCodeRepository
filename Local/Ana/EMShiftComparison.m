niter=20;
nx=100;
sh=10;
prRange=20;
xs=(-nx/2:nx/2-1)';
sigmaG=1;  % SD of the gaussian signal
m0=exp(-xs.^2/(2*sigmaG^2));
data0=circshift(m0,[sh 0]);
emModel=m0;
bayesModel=m0;
nRealiz=1000;
nSet=10;

nDis=100;
sigmaNoise=2;
sigmaModel=sigmaNoise;

emResult=zeros(nx,1);
emSumProb=zeros(nx,nSet);
bayesResult=zeros(nx,1);
bayesSumProb=zeros(nx,nSet);
bayesShifts=zeros(nRealiz,1);

for iRealiz=1:nRealiz
    
    data=repmat(data0,1,nSet)+sigmaNoise*randn(nx,nSet);
    emModel=m0+sigmaModel*randn(nx,1);
    for i=1:niter
        [emModel,emProb]=EM_shift_est(emModel,data,sigmaNoise,prRange);
%         plot(emModel);
%         drawnow;
    end;
%     pause(0.2);
    emReconstruction=emModel;
    bayesModel=emModel;
    bayesReconstruction=zeros(size(bayesModel));
    for i=1:nSet
        [bayesShift,bayesProb(:,i)]=bayes_shift_est(bayesModel,data(:,i),sigmaNoise,prRange);
        bayesShifts(iRealiz,i)=bayesShift;
        bayesReconstruction=bayesReconstruction+circshift(data(:,i),[1-round(bayesShift) 0]);
    end;
    bayesReconstruction=bayesReconstruction/nSet;
    
    emResult=emResult+emReconstruction;
    emSumProb=emSumProb+emProb;
    bayesResult=bayesResult+bayesReconstruction;
    bayesSumProb=bayesSumProb+bayesProb;
    
    if mod(iRealiz,nDis)==0
        subplot(2,1,1);
%         plot([mean(emSumProb,2) bayesSumProb]/iRealiz);
        plot([mean(emSumProb,2)]/iRealiz);
        subplot(2,1,2);
        plot([emResult bayesResult]/iRealiz);
        %         plot(data)
        drawnow;
    end;
end;

subplot(2,1,1);
% plot([emSumProb/nRealiz bayesSumProb/nRealiz]);
plot([mean(emSumProb,2)/nRealiz]);
% legend('EM','Bayes');

subplot(2,1,2);
% plot([emResult/nRealiz bayesResult/nRealiz m0]);
plot([emResult/nRealiz m0]);
% legend('EM','Bayes','Original');


