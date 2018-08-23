function str=ess(n,e)
% Conditional 's' or 'es' character.  
% Avoids text like 'Found 1 particles' by using the construct
% disp(['Found ' num2str(n) ' particle' ess(n)]);
% or sprintf('Found %g particle%s',n,ess(n)) 
% returns 's' if n==1, otherwise ''.
% The second argument, if true, causes 'es' to be generated when n==1.
if nargin<2
    e=false;
end;
if n==1
    str='';
elseif e
    str='es';
else
    str='s';
end;
