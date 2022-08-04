 % FTEx3: complex plot of disc, showing FT shift
%  Click in the left window to move the real-space disc around. q exits.
square=0 ;
ctr=2;

disp('FTEx3: complex plot of FT and shifting');
n=1024;
n1=3*n/8+1;
n3=5*n/8;
ctr=n/2+1; % center point
rd=10;  % radius of disc
rs=3;  % size of square
figure(1);
set(gcf,'menubar','none');
% set(gcf,'doublebuffer','on');
set(gcf,'color',[.3 .3 .3]);
SetComplex;
b=1;
i=0; j=0;
while ~(b==32||b=='q') % space or 'q' exits the loop
    if square
        p=zeros(n);
        p(ctr+i-rs:ctr+i+rs, ctr+j-rs:ctr+j+rs)=1;  % square
    else
%         p=disc(n,rd,[129+i 129+j]); % disc
        p=fuzzymask(n,2,rd,rd/5,[ctr+i ctr+j]); % disc
    end;
    p1=p*255;  % Make the displayed version of p.
    p1(ctr,:)=272;  % Vertical line, plotted as about 0.3 + 0i
    p1(:,ctr)=272; % Horizontal line
    % Show the central 1/4 of p1
    subplot(1,2,1); imac(p1(n1:n3,n1:n3)); axis off;
    title([i j]);
    % Compute the FT of p.
    pf=fftshift(fftn(fftshift(p)));
    subplot(1,2,2);
    imacx(pf,0.4); axis off
    subplot(1,2,1);
    [i,j,b]=Myginput(1);
    i=round(i+n1-1-ctr);
    j=round(j+n1-1-ctr);
end;
% subplot(1,2,2);
% imacs(abs(m));

