function p=PoissonProb(n,mean)

p=mean.^n./factorial(n)*exp(-mean);