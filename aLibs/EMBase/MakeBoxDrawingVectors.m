function [boxesX,boxesY,textX,textY]=MakeBoxDrawingVectors(boxCoords,boxRadius,cornerFraction)
% function [boxesX,boxesY,textX,textY]=MakeBoxDrawingVectors(boxCoords,boxRadius,cornerFraction)
% Create x and y coordinates to draw a set of boxes over an image.
% Creates 10 x nb matrices of coordinates of vectors with NaNs to separate
% them. Also creates locations (center below each box) to draw text.
% boxCoords is an n x 2 array in image coordinates (1-based).
% cornerFraction is =1 to draw a square box, <1 for rounded corners, >1
% for sharp corners.
% Example:
%   [bX,bY,tX,tY]=MakeBoxDrawingVectors(coords,128,0.8); %box size 256, rounded.
%   imags(img);
%   hold on;
%   plot(boxesX,boxesY,'color',[.8 .8 0]); % yellow boxes
%   for i=1:size(textX,1)
%       text(textX,textY,num2str(i));
%   end;
%   hold off;


if nargin<3
    cornerFraction=1;  % simple square
end;

numB=size(boxCoords,1); % Number of boxes to construct
if numB<1
    boxesX=NaN;
    boxesY=NaN;
    return;
end;

% Get the box centers in the display coordinates
bx=boxCoords(:,1);
by=boxCoords(:,2);

% Boxes with modified corners
f=cornerFraction; % <1 for turned-down corners
sqx0=[-f f  1  1 f -f -1 -1 -f NaN]'*boxRadius;
sqy0=[-1 -1 -f f 1 1  f  -f -1 NaN]'*boxRadius;

nL=numel(sqx0); % Total lines to make

boxesX=single(ones(nL,numB))*NaN; % By default, we make everything invisible
boxesY=single(ones(nL,numB))*NaN;

for i=1:numB
        boxesX(:,i)=bx(i)+sqx0; % Center coordinate plus all the box twiddles
        boxesY(:,i)=by(i)+sqy0;
end;
boxesX=boxesX(:);
boxesY=boxesY(:);

if nargout>2 % we're returning text coordinates too.
% The function text() requries doubles for coordinates :(
    textX=double(bx); % location for text below the box
    textY=double(by-boxRadius);
end;
