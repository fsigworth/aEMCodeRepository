% AnaData1
% Get extracted vesicle data from a movie

cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Merged')
m=ReadMRC('sq10_007_Jun25_16.24.33m.mrc');
figure(1);
SetGrayscale;
imacs(m)
cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Micrograph')
sm1=ReadMRC('sq10_007_Jun25_16.24.33ala.mrc');
sm2=ReadMRC('sq10_007_Jun25_16.24.33alb.mrc');
sm1=RemoveOutliers(sm1);
sm2=RemoveOutliers(sm2);

cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/movie_frames/sq10')
mov=ReadMovie('Jun25_16.24.33.tif');

% vesicle 23 coordinates
px=2701;
py=815;
% mark the spot
hold on;
plot(px,py,'y+');
hold off;
drawnow;
%%
n=768;   % extracted image size
n2=n/2;  % half-size of extracted image

mx=ExtractImage(m,[px py]/2+.5,n2);  % extracted merged image
sx1=ExtractImage(sm1,[px py],n);     % extracted low-defocus aligned sum
sx2=ExtractImage(sm2,[px py],n);     % high-defocus aligned sum

% display the sum of the first half of the movie.
mvs=sum(mov(:,:,1:17),3);
imacs(rot90(GaussFilt(RemoveOutliers(mvs),.05),3));
hold on;
plot(px-1,py-65,'y+');  % shift according to padding of micrographs
hold off;
drawnow;

% Create the extracted movie
% The padding includes a shift of 1,65
mvx=zeros(n,n,30,'single');
for i=1:30
    q=rot90(mov(:,:,i),3);
    mvx(:,:,i)=RemoveOutliers(q(px-n2:px+n2-1,py-n2-64:py+n2-65));
end;
imacs(sum(mvx,3));  % show the summed frames

%%
cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1')
mi=ReadMiFile('sq10_007_Jun25_16.24.33mi.txt');
vindex=23;
vesm=meMakeModelVesicles(mi,mi.imageSize,vindex,0,0); % no ctf
c=CTF(mi.imageSize,mi.pixA,mi.ctf(1));  % CTF of the first 17 frames
filtVes=real(ifftn(fftn(vesm).*ifftshift(c)));
ves1=-5*ExtractImage(filtVes,[px py],n);  % factor of 5 seems to scale the model correctly.
imacs(ves1);
save('Vesicle23.mat','mvx','sx1','sx2','mx','ves1');
