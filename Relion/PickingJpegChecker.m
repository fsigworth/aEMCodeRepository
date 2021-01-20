% PickingJpegChecker

jpegDir='Picker_jpegs/';

% load allMis.mat
nmi=numel(allMis);

figure(1);
set(1,'color',[.4 .4 .4]);
set(1,'menu','none');
set(1,'toolbar','none');

ind=1;
ind=MyInput('Starting index ',ind);

while ind<nmi
    mi=allMis{ind};
    allMis{ind}.active=true;
    jpegName=[jpegDir mi.baseFilename '_i0.jpg'];
    if exist(jpegName,'file')
        im=imread(jpegName);
        image(im);
        title(ind);
        fprintf('%d',ind);
    else
        disp([num2str(ind) ' not found: ' jpegName]);
    end
    [x,y,b]=ginput(1);
    switch char(b)
        case 'j' % junk
            allMis{ind}.active=false;
            fprintf(' X\n');
        case 'p' % go to previous
            ind=ind-2;
        case 'q' % quit (and set active=1...)
            break;
    end;
    if allMis{ind}.active
        fprintf('\n');
    end;
    ind=ind+1;
    if ind<1
        beep
        ind=1;
    end;
 end;
            
 disp('Done.');
 
            
