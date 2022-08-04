% FTEx3: complex plot of disc, showing FT shift
%
n=128;
square=0;
c=n/2+1; % center coordinate
cp=[c c]; % center point
rd=5;  % radius of disc
rs=3;  % size of square
figure(1);
% set(gcf,'menubar','none');
set(gcf,'menubar','none','doublebuffer','on','color',[.3 .3 .3]);
% set(gcf,'color',[.3 .3 .3]);
SetComplex;
b=1;
i=0; j=0;
iv=randn; jv=randn;
while ~(b==32|b==121) % space or 'y' exits the loop
    if square
        p=zeros(n);
        p(c+i-rs:c+i+rs, c+j-rs:c+j+rs)=1;  % square
    else
        p=disc(n,rd,[c+i c+j]); % disc
    end;
    p1=p*255;
    p1(c,:)=270;
    p1(:,c)=270;
    subplot(1,2,1); imac(p1); axis off;
    title([i j]);
    %  p=zeros(256); p(129,133)=1;
    pf=fftshift(fftn(fftshift(p)));
    m=pf;
    subplot(1,2,2);
    imacx(m,0.8); axis off;
    drawnow;
%     [ip,jp,b]=ginput(1);
    i=round(i+iv);
    j=round(j+jv);
    iv=0.8*iv+randn(1)-.05*i;
    jv=0.8*jv+randn(1)-.05*j;
    pause(0.2);
end;
% subplot(1,2,2);
% imacs(abs(m));

