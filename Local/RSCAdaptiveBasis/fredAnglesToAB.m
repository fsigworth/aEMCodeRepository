function data_axes = fredAnglesToAB(projDirs) 

%% SET GLOBAL FLAGS
%   DEBUG = 0, silent operation. 
%   DEBUG = 1, various plots during processing.
DEBUG = 0;  

data_axes = [];

% Calculate normals for MATLAB reconstruction
degToRad = pi/180;

eulerAngles = [projDirs(:,1)+90 projDirs(:,2) projDirs(:,3)-90];
eulerAngles = eulerAngles*degToRad;

for i = 1:size(projDirs,1)

  psi = -1*eulerAngles(i,3);
  the = eulerAngles(i,2);
  phi = eulerAngles(i,1);

  T1=[cos(psi)  sin(psi)  0;...
     -sin(psi)  cos(psi)  0;...
      0         0         1];

  T2=[cos(the)	0 -sin(the);...	% from FRED's file
	  0         1 0;...
      sin(the)	0 cos(the)];  

  T3=[cos(phi)  sin(phi) 0;...
     -sin(phi)  cos(phi) 0;...
      0         0        1];
  
  T = T1*T2*T3;  
  T = [T(:,1); T(:,2); T(:,3)];

  if i == 1
    coord_axes = T;
  else
    coord_axes = [coord_axes T]; 
  end

  clear psi phi the
end

%% Create Data Structure for MATLAB reconstruction
data_axes = cat(2,data_axes,coord_axes);

if DEBUG==1
    figure(1000);
    plot_vectors(data_axes);
end

end

function plot_vectors(coord_axes)
vsize=size(coord_axes);
for i=1:vsize(2),
        plot3(coord_axes(7,i),coord_axes(8,i),coord_axes(9,i),'k*');
        hold on;
        plot3([0 coord_axes(7,i)],[0 coord_axes(8,i)],[0 coord_axes(9,i)],'k');
end
plot3(0,0,0,'g+');
pause(0.1)

end

