% PlotComplexExponential


z=0:.01:12;
x=cos(z);
y=sin(z);
plot3(z,x,y,'linewidth',2)
ylabel('Real part');
zlabel('Imaginary part');
xlabel('\theta')
xlabel('theta')
grid on;
hold on;
% plot3(z,x,-1+0*y,'k--')
% plot3(z,1+0*x,y,'k--')
plot3(z,0*x,0*y,'k--');
plot3(0,0,0,'k+','markersize',20);
hold off;
axis equal

%
v=VideoWriter('ComplexExponential1.mov','MPEG-4');
v.Quality=95;
open(v);

for i=0:.4:80
    i1=90-i;
    j1=((i-45)^2-45^2)*.002-.1*i;
    ax.View=[i1,j1];
    frame=getframe(gcf);
    writeVideo(v,frame);
    if i<10
        for j=1:3
    writeVideo(v,frame);
        end;
    end;
end;
close(v);

%%
v=VideoWriter('ComplexExponential2.mov','MPEG-4');
v.Quality=95;
open(v);
for i=0:90
    ax.View=[i1*(90-i)/90 j1*(90-i)/90+i];
% for i=0:90,ax.View=[i1*(90-i)/90)+0*.002*((i-45)^2-45^2) i];
    frame=getframe(gcf);
    writeVideo(v,frame);
end;
close(v);
