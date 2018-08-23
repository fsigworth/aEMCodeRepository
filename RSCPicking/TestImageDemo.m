% TestImageDemo.m

TestFakeMicrograph;

%%
pStack=SimpleExtractor(partImage,mi,48);
vStack=SimpleExtractor(simImage,mi,48);

figure(4);
subplot(1,2,1);
imags(sum(pStack,3));
subplot(1,2,2);
imags(sum(vStack,3));

% Note because of the quirks of ImagicDisplay3 that you have to close these figures if you want to re-use them for
% something else!
figure(5);
ImagicDisplay2(pStack);
figure(6);
ImagicDisplay2(vStack);



