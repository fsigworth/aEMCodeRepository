function path=ViterbiPath(p0,B,A)
% function path=ViterbiPath(p0,B,A)
% given ni states and nt time points, p0 is the ni x 1 initial probability
% vector; B is the ni x nt observable probability, and A is a 1 x nk vector
% representing a row of (if nk=3) the tridiagonal matrix.
% B=[ 1 0 0 0 0 0 0 0
%     1 0 0 0 0 0 0 0
%     0 1 0 0 0 0 0 0
%     0 0 1 0 0 0 0 0
%     0 0 0 1 0 0 0 0
%     0 0 0 1 0 0 0 0
%     0 0 1 0 0 0 0 0
%     0 1 0 0 0 0 0 0 ]';
% B=GaussFilt(B,.2);
% p0=[1 .1 0 0 0 0 0 .1]';
% A=[.3 .4 .3];

[ni,nt]=size(B);
nk=numel(A);
nk1=floor(nk/2); % number of off-diagonal elements on one side
ks=(-nk1:ni-nk1-1)';  % shift of rows

T1=zeros(ni,nt);
T2=ones(ni,nt);


P=p0(:).*B(:,1);
T1(:,1)=P/sum(P);

for i=2:nt
    Q=T1(:,i-1)*A;  % 3-column matrix
    for j=1:nk
        Q(:,j)=circshift(Q(:,j),nk1+1-j);
    end;
    [mxv,k]=max(Q,[],2);  % which column is best?
    P=B(:,i).*mxv;
    T1(:,i)=P/sum(P);
    T2(:,i)=ks+k;
%     plot(T1);
end;
T2=mod(T2-1,ni)+1;  % wrap around, smallest value is 1

[~,p]=max(T1(:,nt));
path=zeros(nt,1);
path(nt)=p;
for i=nt-1:-1:1
    p=T2(p,i+1);
    path(i)=p;
end;
