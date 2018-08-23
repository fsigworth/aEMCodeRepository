function [loc,b]=GetClick(mode,pointerType,ha)
% 
% ha is an Axes handle.
% First, call this to set up the cursor tracking:
%   GetClick('init');
% or for example
%   GetClick('init','circle',gca);
% pointerType can be arrow, circle, cross, crosshair, square.  Default is
% 'cross'.
% Then call this and check for b~=0
%   [loc,b]=GetClick('read');
% Loc is given in the coordinate system of the axes used for 'init'.
% mouse buttons are 1: left button, 2: center, 3: right, 4: double-click
% If a key was pressed, b is class char, containing the key character.
% The returned values are zero unless a button or key press has happened,
% with the cursor within the given axes.



% fs, Sep 2014

persistent pos button key

if nargin<2
    pointerType='cross';
end;
if nargin<3
    ha=gca;
end;


switch lower(mode)
    case 'init'
        if ~isa(ha,'matlab.graphics.axis.Axes')
            warning('h should be an axes handle');
        end;
        hf=ha;
        level=0;
        while ~isa(hf,'matlab.ui.Figure') && level<5 % max. tree size
            level=level+1;
            hf=get(hf,'Parent');
        end;
        ha.ButtonDownFcn=@ButtonDownFcn;
        hf.KeyPressFcn=@KeyFcn;
        hf.WindowButtonMotionFcn=@MotionFcn;
        %         ha.HitTest='on';  % default anyway
        pos=[0 0];
        button=0;
        key=char(0);
        
        %     Create my pointer here
        if strcmpi(pointerType,'square')
            % Pointer map
            val=1;  % 2 to make the box white.
            pn=16; %  Make a 16 x 16 square
            range=1:pn-1;
            ptrMap=NaN+zeros(pn,pn);
            ptrMap(1,range)=val;
            ptrMap(pn-1,range)=val;
            ptrMap(range,1)=val;
            ptrMap(range,pn-1)=val;
            set(hf,'pointer','custom');
            set(hf,'pointer','custom','pointershapecdata',ptrMap,...
                'pointershapehotspot',[1 1]*(pn/2-1));
        else
            set(hf,'Pointer',pointerType);
        end;
        
        
    case 'read'
        %         Check that coordinates are in bounds, otherwise we return
        %         all zeros.
        limits=[ha.XLim' ha.YLim'];
        if any(pos<limits(1,:) | pos>limits(2,:))
            pos=[0 0]; % if not, we zero out all the variables.
            button=0;
            key=char(0);
        end;
        loc=pos;
        if double(key)>0
            b=key;
        else
            b=button;
        end;
        pos=[0 0];
        button=0;
        key=char(0);
end;


function ButtonDownFcn(ha, eventData) % called with an axes handle
%         We update the variables pos and button
pt=ha.CurrentPoint;
pos=pt(1,1:2);
% To find out which button was pressed, have to find the Figure handle
hf=ha;
while ~isa(hf,'matlab.ui.Figure')
    hf=hf.Parent;
end;
%        and then decode the button number.
switch get(hf,'selectiontype')
    case 'normal'
        button=1;
    case 'extend'
        button=2;
    case 'alt'
        button=3;
    case 'open'
        button=4;
    otherwise
        warning('unrecognized selection type');
end;
get(hf,'selectiontype')
end


function KeyFcn(hf,eventData)  % called with a figure handle
key=get(hf,'currentcharacter');  % store the character
if numel(key)<1
    key=char(0);
end;
pt=get(get(hf,'currentaxes'),'currentpoint');
pos=pt(1,1:2);
end

function MotionFcn(hf,eventData)  % called with a figure handle.
get(hf,'currentPoint');  % force the current point to be recorded.
end

end