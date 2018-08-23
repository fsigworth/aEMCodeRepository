function out=Symmetrize3(in,sym)
% function out=Symmetrize3(in,sym)

out=in;
switch sym
    case 2
        out=out+ImgRot90(in,2);
    case 4
        for i=1:sym-1
            out=out+ImgRot90(in,i);
        end;
    otherwise
        error(['Symmetry not supported: ' num2str(sym)]);
end;
