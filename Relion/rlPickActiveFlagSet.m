function j=rlPickActiveFlagSet(si0,txt)
if nargin<2
    txt='Which active flag set';
end;
disp([num2str(size(si0.miIndex,1)) ' particles total.']);
j=size(si0.activeFlags,2);
for i=1:j
    if numel(si0.activeFlagLog)<i
        si0.activeFlagLog{i,1}='';
    end;
    disp([num2str([i sum(si0.activeFlags(:,i))]) '  ' si0.activeFlagLog{i}]);
end;
j=MyInput(txt,j);
end