function [ind,loc,b]=imagsar(imageStack, contrastFraction, showLabels, LabelOffset)
% function imagsar(imageStack, contrastFraction, ShowLabels, LabelOffset)
% Interactive display of an image array.
% --previously ImagicDisplay3--
% Alternative call: recover clicked-on panels by calling with no
%   arguments.
% imagsar('GetClick');
% imagsar('Mark',indices,markInfo);
%  markInfo has the fields marker and markerSize.
% 
% Make an EMAN-like, resizable, scrollable display of the image stack in the current
% figure window (or create one if no figures exist).  Only the first
% argument is required.
% contrastFraction is the fraction of grayscale to skip; .01 means contrast
% is raised such that it saturates with lowest and highest 1% of pixel
% values.  Grayscale scaling is global, not individual for images.
% If ShowLabels=1 only first one in each
% row is labeled.  If ShowLabels=0 (default) no labels are drawn.
% LabelOffset is an offset added to each panel number; by default this is
% zero, so the first panel is panel 1.  Resizing works as
% long as the window exists.

% fs, Sep 2014, based on ImagicDisplay by Yunhui Liu, Aug 2011 based on an earlier version by fs

persistent click pos button col_max n1 n0 marks
q=get(groot,'ScreenSize');
monitorScale=2-(q(3)>1680);  % small monitor means retina display....
% q
% monitorScale

ind=0;
loc=[0 0];
b=0;

if nargin<1
    imageStack='GetClick';
end;
if ischar(imageStack)
    switch lower(imageStack)
        case 'getclick'
            ind=click;
            loc=single(pos);
            b=button;
            button=0;
            click=0;
            pos=[0 0];
            return
        case 'mark'
            inds=contrastFraction;
            markInfo=showLabels;
            for k=1:numel(inds)
                mk=markInfo(min(numel(markInfo),k));
                % Add to the marks list for refreshing the display
%                 km=size(marks,1);
                marks{end+1,1}=inds(k);
                marks{end,2}=mk;
            end;
            DrawMarks(inds,markInfo);
    end;
else  % Initialization: we have a stack to draw
    click=0;
    button=uint8(0);
    pos=single([0 0]);
    col_max=0;
    marks=cell(0,0);
    
    if nargin<4
        LabelOffset=0;
    end;
    if nargin<3
        showLabels=0;    % Labels only at the start of columns
    end;
    if nargin<2
        contrastFraction=0;
    end;
    if contrastFraction >1 || contrastFraction < 0
        contrastFraction
    end;
    [n0, ~ ,nim]=size(imageStack);
    
    [rangeMin,rangeMax]=Percentile(imageStack(:),contrastFraction);
    n1= n0+1;            % add  1 pixel as border
    
    h=gcf;
    hSP=[];
        
    set(h,'Units','pixels','ResizeFcn',@redraw_boxes,'Toolbar','none','ButtonDownFcn',@button_down);
    
    api=struct;
    fsize=9;  % smallest readable size
    
    redraw_boxes(0,0);
    
end; % nargin



    function ButtonDownFcn(src, eventData)
        offs=[3 3];
%         The parent is an axes object.
        pt=get(get(src,'parent'),'currentpoint');
        % To find out which button was pressed, have to find the Figure handle
        h=src;
        while ~isa(h,'matlab.ui.Figure')
            h=get(h,'parent');
            pos=get(h,'children');
        end;
        %        and then decode the button number.
        button=0;
        switch get(h,'selectiontype')
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
%         disp(button)
        pos=pt(1,1:2)-offs;  % update the position
        %         xind=floor((pt(1,1))/n1);
        %         yind=floor((pt(1,2))/n1);
        %         Compute the index of the panel that was clicked on
        ind=floor((pt(1,1))/n1)+col_max*floor((pt(1,2))/n1)+1;
%         if ind>nim
%             return
%         end;
        click=ind;
        %         p=find(clicks==ind,1);
        %         if numel(p)>0 % already selected this one.
        %             clicks(p)=[];
        %             marker='ko';  % draw a black circle
        %             disp(-ind);
        %         else              % a new selection
        %             if numel(clicks)>0
        %                 clicks(end+1)=ind;
        %             else
        %                 clicks=ind;
        %             end;
        %             marker='ro';
        %             disp(ind);    % draw a red circle
        %         end;
        %         plot(xind*n1+n0/2,yind*n1+n0/2,marker,'markersize',n0*mSize);
        %         disp(pos)
    end

    function KeyFcn(src,eventData)
        button=get(src,'currentcharacter');  % store the character
        offs=[3 3];  % coordinate offset
        pt=get(get(src,'currentaxes'),'currentpoint');
        pos=pt(1,1:2)-offs;
        ind=floor((pos(1))/n1)+col_max*floor((pos(2))/n1)+1;
        click=ind;
        %         plot(xind*n1+n0/2,yind*n1+n0/2,marker,'markersize',n0*mSize);
        %         disp([button ' ' num2str(ind)]);
        %         disp(pos)
    end

    function MotionFcn(hFig,eventData)
        get(hFig,'currentPoint');  % force the current point to be recorded.
        %         pta=get(get(hFig,'currentaxes'),'currentpoint')
    end

    function  redraw_boxes(src,evt)
        h0=gcf;
        clf;
        P=get(h0,'Position');
        %     fsize=max(7,min(10,round(showboxsize/4)));  % font size
        win_width = P(3); win_height=P(4);
        col_max= floor(win_width/n1);  % each line can contain as max as col_max number image
        rows = ceil(nim/col_max);
        
        bigmatrix =ones(rows*n1,col_max*n1)*rangeMin; % the whole image
        
        % construct the bigmatrix (assemble the small images as one image)
        index = 1;  % indexes the images
        for j=1:rows
            yOrg1=(j-1)*n1;
            for i= 1:col_max
                xOrg=(i-1)*n1;
                if index<=nim
                    bigmatrix(yOrg1+1:yOrg1+n0,xOrg+1:xOrg+n0) = rot90(imageStack(:,:,index));
                end;
                index =index+1;
            end
        end
        
% hu=uipanel('position',[0 0 .5 1]);
        apos1=[0 win_height-rows*n1 n1*col_max n1*rows];
        ha=axes(h0,'Units','pixels','Position',apos1);  % axes handle
        set(ha, 'color', [0 0 0],'xtick',[],'ytick',[]);
        if max(bigmatrix(:))<=min(bigmatrix(:)) % nothing to show
            return
        end;
        
        hIm = imshow(bigmatrix,'InitialMagnification',100,'DisplayRange',[rangeMin rangeMax]);
        hold on;

        if showLabels>0  % Show some labels
            set(gca,'defaultTextHorizontalAlignment','right',...
                'defaultTextVerticalAlignment','bottom',...
                'defaultTextColor','w','defaultTextBackgroundColor',[.2 0 .2],...
                'defaultTextFontsize',fsize);
            %             ,'defaultTextUnits','pixels');
            index=1;
            for i=1:rows
                text(n1-3 ,i*n1 ,num2str(index+LabelOffset));
                index=index+col_max;
            end
        end
        if (rows*n1 > win_height) % have to make a scroll panel
            hSP = imscrollpanel(h0,hIm);
%             set(hSP,'Units','normalized','Position',[0 0 1 1])
%          apos2=P;
% %         apos2(2)=0;
%          set(hSP,'Units','pixels','Position',apos2)
            api = iptgetapi(hSP);
            api.setMagnification(monitorScale);
            api.setVisibleLocation(0, 0) ;
            api.setImageButtonDownFcn(@ButtonDownFcn);
            set(gcf,'keypressfcn',@KeyFcn);
            set(gcf,'windowbuttonmotionfcn',@MotionFcn);
            % %             api.addNewLocationCallback(@NewLocationFcn);
            %             rect=api.getVisibleImageRect();
            hold on;
        else
            q=get(gca,'children');  % the last of these is the image.
            q=q(end);
            set(q,'buttonDownFcn',@ButtonDownFcn);
%             set(q,'buttonDownFcn',@q.Parent.ButtonDownFcn);
%             set(q(end),'hittest','off');
%             gca.ButtonDownFcn=@ButtonDownFcn;

            %              set(gcf,'buttonDownFcn',@ButtonDownFcn);
        end
        %         api = iptgetapi(hSP);
        %         loc=api.getVisibleLocation()
        DrawMarks(marks);

    end

    function DrawMarks(inds,markInfo)
        if nargin<2 % inds is really the marks cell array
            for jm=1:size(inds,1)
                mk=inds{jm,2};
                ind=inds{jm,1};
                x=mod(ind-1,col_max)*n1+n0/2;
                y=floor((ind-1)/col_max)*n1+n0/2;
                h=plot(x,y,mk.marker,'markersize',mk.markerSize,...
                    'buttondownfcn',@ButtonDownFcn);  % don't capture mouse clicks
            end;
        else
            for jm=1:numel(inds)
                mk=markInfo(min(numel(markInfo),jm));
                x=mod(inds(jm)-1,col_max)*n1+n0/2;
                y=floor((inds(jm)-1)/col_max)*n1+n0/2;
                h=plot(x,y,mk.marker,'markersize',mk.markerSize,...
                    'buttondownfcn',@ButtonDownFcn);  % don't capture mouse clicks
            end;
        end;
    end

end % ImagicDisplay

function [val, upperVal]=Percentile(x,fraction)
%                 val=Percentile(x,fraction)
% [lowerVal upperVal]=Percentile(x,fraction)
% Returns the value of the element of x closest to the given fraction
% (between 0 and 1) of the distribution of x.  X can be of any dimension;
% all elements are used.
% If two output variables are expected, we find the lower and upper values
% corresponding to fraction and (1-fraction).  For example if fraction =
% .01, then the values of the 1% and 99% elements are returned.
x=x(:);
n=numel(x);
xs=sort(x);
nfrac=round(n*fraction)+1;
if nfrac<1 || nfrac>n+1
    error('Fraction out of bounds, must be in [0..1]');
end;
nfrac=min(n,nfrac);  % don't allow it to go beyond n
val=xs(nfrac);

if nargout>1
    upperNFrac=round(n*(1-fraction))+1;
    upperNFrac=min(n,upperNFrac);
    upperVal=xs(upperNFrac);
end;
end 
