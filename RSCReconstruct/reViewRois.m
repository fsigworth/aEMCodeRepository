% reViewRois.m
% Look at roi.mat files and show the class means

nBin=3;
doCTF=1;
k=1e-7;

[names, pa]=uigetfile('*roi.mat','Select roi.mat files','multiselect','on');
if isnumeric(pa)
    return
end;
cd(pa);
if ~iscell(names)
    names={names};
end;

nFiles=numel(names);
fileIndex=1;
b=1;

while b ~= 'q'
    
    if fileIndex>nFiles
        fileIndex=nFiles;
        beep;
    elseif fileIndex<1
        fileIndex=1;
        beep;
    else
        roiName=names{fileIndex};
        
        % Load either roi or gRoi structure.
        disp(['loading ' roiName]);
        roi=load(roiName);
        if isfield(roi,'gRoi')
            roi=roi.gRoi;
        elseif isfield(roi,'roi')
            roi=roi.roi;
            roi
        else
            disp('Wrong structure.');
            roi
            return
        end;
    
        
        % Show the class means
        [nx,ny,nim]=size(roi.classMeans);
        classMeansNorm=zeros(nx,ny,nim,'single');
        if doCTF
            for i=1:nim
                cNorm=ifftshift(roi.classNorms(:,:,i));
                if cNorm(:)'*cNorm(:)>k.^2 %% something there
                    classMeansNorm(:,:,i)=real(ifftn(fftn(roi.classMeans(:,:,i))./(k+cNorm)));
                end;
            end;
        else
            % Simple normalization by std
            cMeans=reshape(roi.classMeans,nx*ny,nim);
            std0=std(cMeans,1);
            sds=max(std0/1000,std(cMeans,1));
            classMeansNorm=reshape(cMeans./repmat(sds,nx*ny,1),nx,ny,nim);
        end;
        
        figure(1);
        ImagicDisplay3(BinImage(classMeansNorm,nBin));
        set(gcf,'name',roiName);
                
    end;  % if fileIndex in range
    
        
    b=0;
    while ~any(b=='npq')
        [partIndex,coords,b]=ImagicDisplay3;
        pause(0.1);
    end;
    
    switch b
        case 'n'
            fileIndex=fileIndex+1;
        case 'p'
            fileIndex=fileIndex-1;
    end;
end;
disp('Done.');


return


%%

% Look at the geometry of YClick vs rVesicle+-mbnOffset
si1=si;
np=numel(si1.miIndex);
rsos=false(np,1);
for i=1:np
    rsos(i)=si1.mi{si1.miIndex(i)}.particle.picks(si1.miParticle(i),7);
end;

sel=~rsos;
mo=si1.mbnOffset;
figure(3);
plot(si1.rVesicle(sel),si1.yClick(sel),'.',si1.rVesicle(sel),si1.rVesicle(sel)+mo,'-');
