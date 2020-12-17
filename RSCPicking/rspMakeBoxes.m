function [boxesX,boxesY,textX,textY]=rspMakeBoxes(mi,dis,boxCoords,cornerFraction)
% dis is the rsp display structure
% boxRadius is half-size of boxes in original pixel size.
% boxCoords is an n x 4 array (x,y,ok,vesicle no.) in original-micrograph coordinates.
% Corner fraction is =1 to draw a square box, 
% If ok<=0 the box at that location isn't plotted.
% Creates 13 x nb matrices of coordinates which create boxes e.g. as
% plot(boxesX,boxesY,'color',[.8 .8 0]);
% Also creates locations (center below each box) to draw text.

if nargin<4
    cornerFraction=1;  % simple square
end;

% rFrac=0.8;
rFrac=0;  % boxes lying farther than rFrac x radius from the center of a
% vesicle get the radial vector shown.

if ~ismatrix(boxCoords)  % handle the case that we've copied from a 3D array
    boxCoords=shiftdim(boxCoords,1);
end;
textX=0;
textY=0;


[numB, ne]=size(boxCoords); % Number of boxes to construct, no. of columns
if numB<1
    boxesX=NaN;
    boxesY=NaN;
    return;
end;
if ne<4
    boxCoords(1,4)=0;  % expand the array with zeros for vesicle indices
end;

% Get the box centers in the display coordinates
bx=boxCoords(:,1)/dis.ds-dis.org(1)+1;
by=dis.size(2)-boxCoords(:,2)/dis.ds+dis.org(2);

% rb=boxRadius/dis.ds;  % 'radius' of a box

rb=ceil(dis.currentBoxSize/(dis.ds*mi.pixA*2));  % 'radius' of a box

textX=double((bx-rb).*(boxCoords(:,3)>0)); % location for text below the box
textY=double((by-rb).*(boxCoords(:,3)>0));

% Boxes with modified corners
f=cornerFraction; % <1 for turned-down corners
sqx0=[-f f  1  1 f -f -1 -1 -f NaN 0 0 NaN]'*rb;
sqy0=[-1 -1 -f f 1 1  f  -f -1 NaN 0 0 NaN]'*rb;

pv1=numel(sqx0)-2; % pointer to where we'll fill in coordinates of the
% vector that we draw to the vesicle center
pv2=pv1+1;

% 
% fr=0.8;
% sqx0=[-rb*fr rb*fr    rb   rb   -rb*fr -rb NaN 0 0 NaN]';
% sqy0=[-rb     -rb  -rb*fr rb*fr    rb  rb*fr -rb NaN 0 0 NaN]';

nL=numel(sqx0); % Total lines to make

boxesX=single(ones(nL,numB))*NaN; % By default, we make everything invisible
boxesY=single(ones(nL,numB))*NaN;

for i=1:numB
    if boxCoords(i,3)>0 % particle is ok
        boxesX(:,i)=bx(i)+sqx0; % Center coordinate plus all the box twiddles
        boxesY(:,i)=by(i)+sqy0;
        vIndex=boxCoords(i,4);
        if vIndex>numel(mi.vesicle.x)
            warning('Inconsistent vesicles and picks, can''t draw boxes');
            return;
        end;
        if vIndex>0 % a vesicle is marked, draw the vector toward the center
            % Get the displacement from box center to vesicle center
            vx=mi.vesicle.x(vIndex)/dis.ds-dis.org(1)+1-bx(i);
            vy=dis.size(2)-mi.vesicle.y(vIndex)/dis.ds+dis.org(2)-by(i);
            dist=hypot(vx,vy); % length of the full vector
            if dist>rFrac*mi.vesicle.r(vIndex)/dis.ds % if the vector is long enough, draw it.
                rbe=rb/max(abs(vx),abs(vy));
                s1=rbe; % Vector start point
                s2=1;  % Vector end: all the way to the center
                boxesX(pv1:pv2,i)=[s1 s2]'*vx+bx(i);
                boxesY(pv1:pv2,i)=[s1 s2]'*vy+by(i);
            end;
        end;
    end;
end;
boxesX=boxesX(:);
boxesY=boxesY(:);
