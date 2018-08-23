


function [P_ws,dW2,Flag]=P_wc2(alphas,betas,d_a,d_b,sigmaB,sigmaR,m,r,y,sigmaC)
% Make a map for p(w,beta,alpha,sigmaB,SigmaR) 
% Input: vector,blobing, rocking and pixel size of particle image


% for test: 
% alphas=0;
% beta=pi*45/180;
% sigmaB=2;
% sigmaR=1;
% p=P_ab(alpha,beta,sigmaB,sigmaR);
% or
% alphas=(-10:5:10);
% betas=(0:10:180)/180*pi; 

% make a sphere map for alpha=0 
%figure(12);
%SetGrayscale;



sigmaC2=sigmaC^2;
%r1=128; % radius
%r=r1;
r1=r;
%d=r1*4; % size of matrix
%pc=GaussMask(2*m,2,sigmaC);

nalpha=numel(alphas);
nbeta=numel(betas);
ngamma=9;
Flag=zeros(nalpha,nbeta); % whether need check certain w
dW2=Flag;
%p_w=zeros(d,2);
P_ws=zeros(m,m,nalpha,nbeta);

for i=1:nbeta
    beta=betas(i);
    r_beta=r1*sin(beta);
% define the center and std for the gauss distribution
    %ctr=[d/2+1,r_beta+d/2+1];
    if abs(beta-pi/2)<8*pi/18
        sigmaX2=d_a^2*r_beta^2+sigmaR^2*r_beta^2; 
        sigmaY2=d_b^2*(cos(beta))^2*r1^2+sigmaR^2*(cos(beta))^2*r1^2+sigmaB^2*(sin(beta))^2; % for testing
        %sigmaW=sqrt([sigmaX2,sigmaY2]);

        for j=1:numel(alphas)
            alpha=alphas(j);

            rotM=[cos(alpha),-sin(alpha);sin(alpha),cos(alpha)]; % counter clockwise rotate by alpha
            %SigmaC=rotM*Sigma;
            dxC=r1*sin(beta)*sin(-alpha)-0;  % the origin was counter clockwise by -alpha
            dyC=r1*sin(beta)*cos(-alpha)-y;
            dCT=rotM*[dxC dyC]'+m+1; % coordinate from the center, rotate back, counter clockwise  by alpha
           % rot2M=[cos(alpha),-sin(alpha);sin(alpha),con(alpha)];
           % dCT2=rot

            dx=(dCT(1)*sigmaC2+(m+1)*sigmaX2)/(sigmaC2+sigmaX2);
            dy=(dCT(2)*sigmaC2+(m+1)*sigmaY2)/(sigmaC2+sigmaY2);
            sigma_X2=sigmaC2*sigmaX2/(sigmaX2+sigmaC2);
            sigma_Y2=sigmaC2*sigmaY2/(sigmaY2+sigmaC2);
            Residual=-dx^2/sigma_X2-dy^2/sigma_Y2+dCT(1)^2/sigmaX2+dCT(2)^2/sigmaY2+2*(m+1)^2/sigmaC2;
            if (-Residual)<-25% Value too small
            else
                Flag(j,i)=1;
                p_w=exp(-1/2*Residual)*GaussMask(2*m,2,sqrt([sigma_X2,sigma_Y2]),[dx,dy])/(4*pi^2*sigmaC2*(sigma_X2+sigma_Y2));
                P_ws(:,:,j,i)=p_w(m/2+1:3*m/2,m/2+1:3*m/2);
                dW2(j,i)=dxC^2+dyC^2;
            end;
                
        end;
    else
        std=r1*pi/12;
        if abs(y)<m
            pc=GaussMask(m,2,sigmaC)/(2*pi*sigmaC2);
            p_w=GaussMask(2*m,2,std,[0,-y]+m+1)/(2*pi*std^2);
            P_ws(:,:,1,i)=p_w(m/2+1:3*m/2,m/2+1:3*m/2).*pc;
            dW2(1,i)=y^2;
            Flag(:,i)=1;
        end;
    end;
%      subplot(4,5,i);
%      imacs(P_ws(:,:,1,i));
% 
%      title(num2str(beta/pi*180));
%      drawnow;
end;
   % pbw=ImageArray(P_ws(:,:,1,:));
   % imacs(pbw);
%drawnow; 
%rotate the mask to the desired alpha
Pabr=1/((nalpha*(nbeta-2)+1)*ngamma);
P_ws=Pabr*P_ws;
end