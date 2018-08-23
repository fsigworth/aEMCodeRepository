function imovie(m, pauseTime, repeats, globalScale)
% function imovie(m, pauseTime, repeats)
% Show an nx x ny x nim array as a movie of nim frames
% the array may be >3 dimensional.  pauseTime is the time that each frame
% is displayed.
if nargin<2
    pauseTime=.01;
end;
if nargin<3
    repeats=1;
end;
if nargin<4
    globalScale=1;
end;

[nx, ny, nim]=size(m);

if globalScale
    m=imscale(m,256,.0001);
end;    
for j=1:repeats
for i=1:nim
    if globalScale
        imaga(m(:,:,i));
    else
        imags(m(:,:,i));
    end;
    title(i);
    pause(pauseTime);
end;
end;