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
axis equal
ax.View=[90 0];
drawnow;
for i=1:90
    ax.View=[90-i,((i-45)^2-45^2)*.002];
    axis equal
    drawnow;
    pause(.1)
end;
%%
ax.View=[45 0];
%%
for i=0:90,ax.View=[.004*((i-45)^2-45^2) i];
    drawnow;
    
pause(0.05)
end;
