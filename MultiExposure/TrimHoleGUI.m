function [] = GUI_13()
% Demonstrate how to display & change a slider's position with an edit box.  
% Slide the slider and it's position will be shown in the editbox.  
% Enter a valid number in the editbox and the slider will be moved to that 
% position.  If the number entered is outside the range of the slider, 
% the number will be reset.
%
%
% Author:  Matt Fig
% Date:  7/15/2009
ix=100;
iy=500;
iw=1300;
ih=600;
S.fh = figure('units','pixels',...
              'position',[ix iy iw ih],...
              'menubar','none',...
              'name','GUI_13',...
              'numbertitle','off',...
              'resize','off');
S.ax = axes('units','pixels','position',[40 40 iw-200,ih-50]);
plot(S.ax,sin(0:.1:10));

          for i=1:3
    S.sl(i) = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[iw-120 10+50*(i-1) 100 30],...
                 'min',0,'max',100,'val',0);             
end;
set(S.sl,'call',{@ed_call,S});  % Shared Callback.


function [] = ed_call(varargin)
% Callback for the edit box and slider.
[h,S] = varargin{[1,3]};  % Get calling handle and structure.

for i=1:3
    v(i)=get(S.sl(i),'value');
end;
x=0:.001:1;
plot(S.ax,sin(x*v(1)/2+x.^2*v(2)));

% S
% h
% switch h  % Who called?
%     case S.ed
%         S.ed
%         L = get(S.sl(1),{'min','max','value'});  % Get the slider's info.
%         E = str2double(get(h,'string'));  % Numerical edit string.
%         if E >= L{1} && E <= L{2}
%             set(S.sl(1),'value',E)  % E falls within range of slider.
%         else
%             set(h,'string',L{3}) % User tried to set slider out of range. 
%         end
%     case S.sl
%         set(S.ed,'string',get(h,'value')) % Set edit to current slider.
%     otherwise
%         % Do nothing, or whatever.
% end