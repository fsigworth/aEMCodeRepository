% Element-wise l_1 or l_0 shrinkage on cell
function coefs = prox_regularizer(coefs, lambda, hard)
if(iscell(coefs))
  if(hard)
    coefs = cellfun(@(f)wthresh(f,'h',lambda), coefs, 'Un', 0);
  else
    coefs = cellfun(@(f)wthresh(f,'s',lambda), coefs, 'Un',0);
  end
else
  if(hard)
    coefs = wthresh(coefs,'h',lambda);
  else
    coefs = wthresh(coefs,'s',lambda);
  end
end
