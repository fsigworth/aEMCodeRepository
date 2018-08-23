% CCX3.m running the exchange
N=2e6;
doPrint=0;
showAllAlphaBeta=0;

aC=AlphaLabs;
bC=BetaLabs;
pC=Carpathia;


% We should automatically discover fields, but here it's hardwired.
allVars=zeros(3,N);
allOks=false(3,N);

maxVals=[5 5 2.5]'; % log maximum values
for i=1:N
     vars=round(10.^(maxVals.*rand(3,1)));
%      vars=round(1e6*rand(3,1));
    aC.aPrice=vars(1);
    pC.aPrice=aC.aPrice;
    
    aC.bPrice=vars(2);
    bC.bPrice=aC.bPrice;
    pC.bPrice=aC.bPrice;
    
    aC.volume=vars(3);
    bC.volume=aC.volume;
    pC.volume=aC.volume;
    
%     aC.altCompensation=vars(4); % leave these at zero
%     pC.altCompensation=vars(4);
    
    aOk=AlphaLabs(aC);
    bOk=BetaLabs(bC);
    pOk=Carpathia(pC) || showAllAlphaBeta;
% pOk=true;
    allVars(:,i)=vars;
    allOks(:,i)=[aOk bOk pOk]';
    
    if aOk && bOk && pOk && doPrint
        disp(vars');
    end;
end;    
%%
oks=all(allOks(1:3,:),1);
nOks=sum(oks);
okPars=allVars(:,oks);

[sortVols,sortInds]=sort(okPars(3,:),'descend');

nPick=ceil(min(nOks,max(10,nOks/10))); % pick top 10%
pickPars=okPars(:,sortInds(1:nPick)); 

centroid=mean(pickPars,2);
dists=sqrt((pickPars(1,:)-centroid(1)).^2 + (pickPars(2,:)-centroid(2)).^2);
[minDist,ind]=min(dists);
finalPars=pickPars(:,ind);

%%
% oks=squeeze(allOks(3,:));

figure(1);
clf;
subplot(122);
plot3(allVars(1,oks),allVars(2,oks),allVars(3,oks),'bo');
xlabel('Asterix price');
ylabel('Betafix price');
zlabel('volume');
grid on
hold on;
plot3(finalPars(1),finalPars(2),finalPars(3),'r.','markerSize',25);
hold off;

subplot(121);
abOks=all(allOks(1:2,:),1);
plot3(allVars(1,abOks),allVars(2,abOks),allVars(3,abOks),'bo');
xlabel('Asterix price');
ylabel('Betafix price');
zlabel('volume');
% plot3(finalPars(1),finalPars(2),finalPars(3),'k.','markerSize',100);
grid on
% hold on;
% % plot3(allVars(1,abOks),allVars(2,abOks),allVars(3,abOks),'bo');
% % xlabel('Asterix price');
% % ylabel('Betafix price');
% % zlabel('volume');
% plot3(finalPars(1),finalPars(2),finalPars(3),'k.','markerSize',100);
% hold off;



%
disp(' ');
disp(' ');
euro=char(8364);
disp('---Alpha contract---')
disp(['Asterix price ' euro num2str(finalPars(1))]);
disp(['Volume ' num2str(finalPars(3))]);
disp('Supply chain violation penalty = ?300000');
disp('Payer provides Zurich warehouse = true');
disp(' ');
disp('---Beta contract---')
disp(['Betafix price ' euro num2str(finalPars(2))]);
disp(['Volume ' num2str(finalPars(3))]);
disp(' ');
disp('---Carpathia contract---');
disp(['Total price ' euro num2str(ceil(sum(finalPars(1:2))*1.01))]);
disp(['Volume ' num2str(finalPars(3))]);
disp(['Supply chain violation penalty = ' euro '300000']);
disp('Payer provides Zurich warehouse = true');
disp(' ');
disp('---CCX revenue--');
cut=(ceil(sum(finalPars(1:2))*1.01))-sum(finalPars(1:2));
disp([euro num2str(cut) ' per unit']);
disp(' ');







