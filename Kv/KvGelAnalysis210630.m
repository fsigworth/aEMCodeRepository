% KvGelAnalysis

% analyze the gel image from the email of Yangyang, 
imc=IMG_5782Rot(2201:2800,601:1900,:);
img=single(imc(:,:,1:2));
img=sum(img,3);
bands=[89 190; 311 430; 563 679; 829 926 ; 889 917; 1090 1172];
nb=size(bands,1);

lanes=zeros(600,6);
for i=1:nb
q=mean(img(:,bands(i,1):bands(i,2)),2);
lanes(:,i)=q;
end;

mysubplot(211);
imags(img)

mysubplot(212)

plot(-lanes(:,1:4))
hold on;
plot(-lanes(:,6)*.1-280,'k-');
legend('Kv','7 days','lipo','pellet','stds');

print('fig.jpg','-r300','-djpeg');

