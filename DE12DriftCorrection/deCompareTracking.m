% deCompareTracking
tfs=14;
doCrop=1;
% doCrop=0;

cd('/Volumes/TetraData/EMWork/Anchi/12apr26c');
baseFilename='12apr26c_b_00029gr_00011sq_v02_00006hl_v01_00008ed';
if doCrop 
    postscr='g2crop-res';
    postscr='g1cropwt05-res';
else
    postscr='g2-res';
    postscr='g2-res';
end;

label=[baseFilename(1:8) baseFilename(46:50) postscr(1:2)];

if doCrop
    xLabel=[label ' (crop) X-shift'];
    xOffset1=.3;
    xOffset2=1;
else
    xLabel=[label '  X-shift'];
    xOffset1=0;
    xOffset2=0;
end;

shx=[];
shy=[];
for i=1:2
    %     load(['Tracking/' baseFilename num2str(i) 'g2-res.mat']);
    load(['Tracking/' baseFilename num2str(i) postscr '.mat']);
    if i==2 && nsh>numel(res.shiftX)
        shx(nsh,:)=[];
        shy(nsh,:)=[];
    end;
    nsh=numel(res.shiftX);
    shx(:,i)=res.shiftX;
    shy(:,i)=res.shiftY;
end;

sigmaX=std(shx(2:nsh,1)-shx(2:nsh,2))
sigmaY=std(shy(2:nsh,1)-shy(2:nsh,2))


nf=nsh;
x=zeros(nf,2);
x(:,1)=(1:nf)';
x(:,2)=xOffset1+x(:,1);
subplot(1,2,1);
plot(x(:,1),shx(:,1),'.-',x(:,2),shx(:,2),'.-','markersize',14);
% ,x(:,1),shx(:,2)-shx(:,1),'.-',...
%      x(:,1),0*x(:,1),'k-');
title([xLabel ' SD= ' num2str(sigmaX)],'fontsize',tfs);
legend({'Even frames','Odd frames'});
xlabel('Averaged frame');
ylabel('Shift, pixels');

subplot(1,2,2);
x(:,2)=xOffset2+x(:,1);
plot(x(:,1),shy(:,1),'.-',x(:,2),shy(:,2),'.-','markersize',14);
% ,x(:,1),shy(:,2)-shy(:,1),'.-',...
%     x(:,1),0*x(:,1),'k-');
title(['Y-shift SD= ' num2str(sigmaY)],'fontsize',tfs);
xlabel('Averaged frame');
ylabel('Shift, pixels');

%%
outName=['Jpeg/' baseFilename postscr '-compare.jpg'];
set(gcf,'paperpositionmode','auto');
print('-djpeg','-r300',outName);
disp(outName);


%         comparison
% g2crop: frameDose=0.99

% sx=.34 when divided by sqrt(2)
% sy=.27 "
%
% g2 no cropping
% sx=.146 when divided by sqrt(2)
% sy=.19  %%this becomes .113 with the new GetWeighting filter with fc=.05
% and exponent is .4 not .1.
