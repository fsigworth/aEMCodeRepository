% CCX2.m running the exchange
N=1e7;
maxTestPrice=8000;  % max price for Asterix or Betafix

% Initialize the structures
aC=Asterix2;
bC=Betafix2;

% pick up the labels
pC=struct;
pC.aPrice=aC.aPrice;
pC.bPrice=bC.bPrice;
pC.volume=aC.volume;

% We should automatically discover fields, but here it's hardwired.
allOks=false(N,3);

maxLogPrice=log10(maxTestPrice);
maxVals=repmat([maxLogPrice maxLogPrice 6],N,1); % log maximum values
    vars=round(10.^(maxVals.*rand(N,3))); % integers
    
    vars(:,3)=5000;  % fix the volume here.
    
    aC.aPrice.val=vars(:,1);
    pC.aPrice.val=vars(:,1);
    
    aC.bPrice.val=vars(:,2);
    bC.bPrice.val=vars(:,2);
    pC.bPrice.val=vars(:,2);
    
    aC.volume.val=vars(:,3);
    bC.volume.val=vars(:,3);
    pC.volume.val=vars(:,3);
    
%     aC.altCompensation=vars(4); % leave these at zero
%     pC.altCompensation=vars(4);
    
    aOk=Asterix2(aC);
    bOk=Betafix2(bC);
    abOk=aOk & bOk;
    
    pC.aPrice.val=vars(abOk,1);
    pC.bPrice.val=vars(abOk,2);
    pC.volume.val=vars(abOk,3);
    [minPrice,ind]=Payer2(pC);
 
    disp(['Best price = ' num2str(minPrice)]);
    disp([pC.aPrice.string ' = ' num2str(pC.aPrice.val(ind))]);
    disp([pC.bPrice.string ' = ' num2str(pC.bPrice.val(ind))]);
    disp([pC.volume.string ' = ' num2str(pC.volume.val(ind))]);
%     if aOk && bOk && pOk && doPrint
%         disp(vars');
%     end;
%     allVars=vars;

%%
% oks=all(allOks(1:3,:),1);
% oks=squeeze(allOks(3,:));
figure(1);
clf;
plot(pC.aPrice.val,pC.bPrice.val,'bo');
hold on;
plot(pC.aPrice.val(ind),pC.bPrice.val(ind),'ro');
xlabel('Asterix price');
ylabel('Betafix price');
hold off;

% subplot(2,2,1);
% plot(allVars(1,oks),allVars(2,oks),'bo');
% xlabel('Asterix price');
% ylabel('Betafix price');
% 
% subplot(2,2,2);
% plot(allVars(1,oks),allVars(3,oks),'bo');
% xlabel('Asterix price');
% ylabel('Volume');
% 
% subplot(2,2,3);
% plot(allVars(2,oks),allVars(3,oks),'bo');
% xlabel('Betafix price');
% ylabel('Volume');
% 

