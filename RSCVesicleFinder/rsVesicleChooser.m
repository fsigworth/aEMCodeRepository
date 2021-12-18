% rsVesicleChooser.m
% Read an mi file, load an image containing refined vesicles.
% Select a subset of vesicles based on clicks.
% Click to select; click again to deselect.
% q: quit and save
% z: quit and don't save.

ds=4;

disp('Getting mi.txt files');
[nm,pa]=uigetfile('*mi.txt','Get mi files','multiselect','on');
if isnumeric(pa) % user clicked Cancel
    return
end;
[basePath,infoPath]=ParsePath(pa);
cd(basePath);
if ischar(nm)
    nm={nm};
end;
%%


for fileInd=1:numel(nm)

    mi=ReadMiFile([infoPath nm{fileInd}]);
    [m,M4,ok,rawImg]=meLoadNormalizedImage(mi,mi.padImageSize/ds,'m');
    % origMicrographXY=M*[xOut;yOut;1];
    % [xOut;yOut;1]=M\[origX origY 1];
    %  assuming zero-based coordinates in each case.
    nv=numel(mi.vesicle.x);
    globalXY=[mi.vesicle.x mi.vesicle.y ones(nv,1,'single')]';
    vesXY=double(1+M4\globalXY); % convert to window coordinates
    allInds=find(mi.vesicle.ok(:,1));
    figure(1); clf;
    imags(GaussFilt(m,.1));
    hold on
    % Show all the vesicles
    ShowCircles(mi,M4,allInds,vesXY,[.8 .5 0]);
    inds=find(mi.vesicle.ok(:,4))';
    % Show our selected ones
    ShowCircles(mi,M4,inds,vesXY,[.2 .5 1]);
    hold off;
    title(nm,'interpreter','none');

    b=0;
    while (b~='q' && b~='z')
        [ix,iy,b]=ginput(1);
        % Find the nearest vesicle
        dists=sqrt(sum([ix-vesXY(1,:)' iy-vesXY(2,:)'].^2,2));
        [minDist,ind]=min(dists);
        p=(inds==ind);
        if any(p)
            disp('delete');
            inds(p)=[];
        else
            inds=[inds ind]
        end;

        imags(GaussFilt(m,.1));

        hold on
        % Show all the vesicles
        ShowCircles(mi,M4,allInds,vesXY,[.8 .5 0]);
        % Show our selected ones
        ShowCircles(mi,M4,inds,vesXY,[.2 .5 1]);
        % Show the indices
        ni=numel(inds);
        for i=1:ni
            j=inds(i);
            text(vesXY(1,j),vesXY(2,j),num2str(j),'fontweight','normal',...
                'color','y','fontsize',12);
        end;
        hold off
        title([nm ':  ' num2str(ni) ' selected'],'interpreter','none');

    end;

    if b=='q'
        mi.vesicle.ok(:,4)=false;
        mi.vesicle.ok(inds,4)=true;
        name=WriteMiFile(mi);
        disp([name ' written']);
        title([nm ' saved.']);

    else
        disp('mi file not saved.')
    end;
end; % for fileInd
return


function ShowCircles(mi,M,inds,vesXY,color)
if nargin<5
    color=[.2 .5 1];
end;
for i=1:numel(inds)

    r1=mi.vesicle.r(inds(i),:)/M(1,1);
    if r1(1)<1
        continue
    end;
    [x,y]=CircleLineSegments(r1,min(10,100/r1(1)));
    x=x+vesXY(1,inds(i));
    y=y+vesXY(2,inds(i));
    plot(x,y,'color',color);
    %                 elseif badVes(i)
    %                     plot(x,y,'r-','HitTest','off');
    %                 end;
end;
end
