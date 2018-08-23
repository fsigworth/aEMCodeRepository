function Render3D(map, threshold, color)
% function Render3D(map, threshold,color)
% Make a 3D rendering of the given map.  The default color is blue.

if nargin<3
    color='blue';
end;
% 'patch'
p = patch(isosurface(map, threshold), 'FaceColor', color, 'EdgeColor', 'none');
% 'isonormals'
isonormals(map,p)
% 'view'
view(3); daspect([1 1 1]); axis equal
% 'camlight'
camlight;  camlight(-80,-10); lighting phong;
% 'done'
drawnow;