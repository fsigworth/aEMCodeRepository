function [x,y,b] = MyBusyread(mode)
% MyBusyread('init'); Call this first, with the relevant window and axes
% selected.
% [x,y,b]=MyBusyread % gets x,y from the current axes and b from the
% current figure. b=1,3 for mouse buttons; ascii code for character.
% arrow keys appear as double(b)=28-31 (left, right, up, down)
% control keys are not interpreted, e.g. ctrl-C gives double(b)=3.
%
% MyBusyread('stop') unhooks the callback function
% fs. Based on imagsar's interactive interface.


persistent pos button

if nargin<1
    mode='GetClick';
end;
switch lower(mode)
    case 'getclick'
        %  disp([pos button]);
        x=single(pos(1));
        y=single(pos(2));
        b=button;
        button=0;
        pos=[0 0];
        
    case 'init'
        button=uint8(0);
        pos=single([0 0]);
        h=gcf;
        set(h,'keypressfcn',@KeyFcn);
        
        set(gca,'ButtonDownFcn',@button_down);
        set(gcf,'WindowButtonMotionFcn',@MotionFcn);
    case 'stop'
        set(gca,'ButtonDownFcn','');
        set(gcf,'WindowbuttonMotionFcn','');
    otherwise
        warning(['Unrecognized argument: ' mode]);
end;


    function button_down(src, eventData)
%         offs=[3 3];
        offs=[0 0];
        pt=get(gca,'currentpoint');
        button=0;
        switch get(gcf,'selectiontype')
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
        pos=pt(1,1:2)-offs;  % update the position
    end

    function KeyFcn(src,eventData)
        button=get(src,'currentcharacter');  % store the character
%         offs=[3 3];  % coordinate offset
        offs=[0 0];  % coordinate offset
        pt=get(get(src,'currentaxes'),'currentpoint');
        pos=pt(1,1:2)-offs;
    end

    function MotionFcn(hFig,eventData)
        get(hFig,'currentPoint');  % force the current point to be recorded.
    end


end

