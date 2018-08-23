function H=rlUnblurFilter(movieSize,frameDose,pixA,kV)
% function H=rlUnblurFilter(movieSize,pixA,kV)
% Obtain the "noise transfer function" of the effective filtering as
% implemented in Grant & Grigorieff's Unblur program.
% frameDose is a scalar giving e/A^2 per frame.
% We evaluate
% N_c ? ak^b + c, with a = 0.245, b = ?1.665, and c = 2.81 at 300kV
% q(k,N)=exp(-N/(2N_c(k)))
% H(k)=sum_i {q(k,Ni)} / sqrt(sum_i {q^2(k,Ni)})
if nargin<4
    kV=300;
end;

n=movieSize(1:2);
nim=movieSize(3);
% doses=cumsum([0:nim-1]'*frameDose);  % dose at the beginning of each exposure

freqs=RadiusNorm(n)/pixA;
critDose=(.245*freqs.^-1.665+2.81)*kV/300; % critical dose as a function of 2D frequency
q1=exp(-frameDose./(2*critDose)); % q for a single exposure
denom=sqrt((1-q1.^(2*nim))./(1-q1.^2)); % =sqrt( 1 + q1^2 + q1^4... + q1^(2[nim-1]) )
num=(1-q1.^nim)./(1-q1);
H=num./denom;
% H=critDose;