function newStack=spTransformStack(stack,sp)
n=size(stack,1);  % get the stack dimension
ds=sp.boxSize/n;  % stack downsampling ratio
npar=size(sp.trans,1);
xytfr=[sp.trans/ds sp.rot single(sp.flip)...
    single(sp.class) zeros(npar,1)];
newStack=TransformImages(stack, xytfr);
