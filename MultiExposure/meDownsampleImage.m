function [mOut,M]=meDownsampleImage(m1,M1,nOut);
% Downsample the image m1 to make an image of size nOut (2x1 vector) while
% handling the scale/shift matrices.
% If nOut doesn't match the aspect ratio, the relatively larger dimension
% is used to set the downsampling factor.
% M1 is the affine transform m1->original micrograph of size mi.imageSize.
% Leave it empty if there is none.
% M is the affine transform mOut->original micrograph
    nIn=size(m1);
    ds1=M1(1,1);
    ds2=max(nIn./nOut); % we pick the dimension best matched by the target.
    mOut=DownsampleGeneral(m1,nOut,1/ds2);
    shift=floor(nIn./ds2-nOut)/2; % The shift from the crop operation
    M2=[ds2 0 -shift(1); 0 ds2 -shift(2); 0 0 1];
    if numel(M1)>0
        M=M1*M2; % Composite matrix mapping output to original image
    else M=M2;
    end;