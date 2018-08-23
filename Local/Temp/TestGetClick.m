% TestGetClick

figure(1);
subplot(2,1,1);
plot(randn(128,1),'hittest','off');
subplot(2,1,2);
h=imacs(randn(128));
set(h,'hittest','off');
GetClick('init');

b='a';

while b~='q'
    [pos,b]=GetClick('read');    
    while b==0
        [pos,b]=GetClick('read');
        pause(0.1);
    end;
    disp([num2str(pos) '  ' num2str(b)]);
end;
