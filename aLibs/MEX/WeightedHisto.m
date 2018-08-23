function [hist, norm]=WeightedHisto(indices, weights, nbins)
% function [hist norm]=WeightedHisto(indices, weights, nbins)
% Uses all elements of 'indices' to direct the summing of all elements of
% all dimensions of 'weights' into bins.
% The returned 'norm' array gives the number of entries summed into
% each bin, which can be used for normalization. 'norm' is therefore a
% histogram of the 'indices'.  Indices out of bounds (<1, >nbins) are not used
% at all.
% This function calls the lower-level mex function WtHist.

ind=int32(indices(:)); % bounds checking is now done in mex file.
nb=int32(nbins);
% if numel(weights) ~= numel(ind)
%     error('indices and weights must have the same number of elements.');
% end;
if isreal(weights)
    [hist, inorm]=WtHist(ind, double(weights(:)), nb);
    if nargout>1
        norm=double(inorm);
    end;
else  % m-code for WtHist
    hist=zeros(nbins,1);
    norm=zeros(nbins,1);
    if nbins < 1
        return
    end;
    
    for i=1:numel(ind)
        k=min(nbins,max(1,ind(i)));
        hist(k)=hist(k)+weights(i);
        norm(k)=norm(k)+1;
    end;
end;

