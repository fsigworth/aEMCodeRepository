% rlExtractClassStack.m
% used to make particle stack from one class, now stored on
% drobo4/relion_sim/

% given a particles.star file, e.g. from a selection step, put all the 
% particle images into a stack.

% classes=13;
inPath='Select/job148/';
classAvg=ReadMRC([inPath 'class_averages.mrcs'],1,1);
[nm,da]=ReadStarFile([inPath 'particles.star']);
%%
pixA=da{1}.rlnImagePixelSize(1);
d=da{2};
np=numel(d.rlnCoordinateX);
np=10000
disp([num2str(np) ' particles']);
for i=1:np
    [slice,fName]=rlDecodeImageName(d.rlnImageName{i});
    m=ReadMRC(fName,slice,1);
    if i==1
        n=size(m);
        stack=zeros([n np],'single');
        ctfs=stack;
    end;
    if mod(i,1000)==0
        fprintf('.');
    end;
    stack(:,:,i)=circshift(m,round([d.rlnOriginXAngst(i) d.rlnOriginYAngst(i)]/pixA));
    pars=rlCTFParsFromStar3Line(nm,da,i);   
    ctfs(:,:,i)=CTF(n,pixA,pars);
end;
fprintf('\n');

disp('Rotating...');
rotStack=rsRotateImage(stack,-d.rlnAnglePsi(1:np));
disp(' done.');
imags(sum(rotStack,3));

save('particles.mat','classAvg', 'stack', 'rotStack', 'ctfs', 'pixA','-v7.3');
