function [loc,b]=WaitForClick
% Like GetClick('read') except it waits for a click before returning.

    [loc,b]=GetClick('read');    
    while b==0
        [loc,b]=GetClick('read');
        pause(0.1);
    end;


