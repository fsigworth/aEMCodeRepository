

%
%     disp('Final subtraction of vesicles');
%     vm=meMakeModelVesicles(mi,n);
%     mv=m-vm;
%     %%
%     figure(2);
%     imacs(GaussFilt(mv,.2));
%
%     figure(1);
%     subplot(2,3,1);
%     imacs(BinImage(m,4));
%     subplot(2,3,2);
%     imacs(vm);
%     subplot(2,3,3);
%     plot(mi.vesicle.r*mi.pixA,mi.vesicle.s,'k.','markersize',10);
%
%     subplot(2,3,4);
%     imacs(BinImage(mv,4));
%     title('Subtracted');
%     subplot(2,3,5);
%     hist(mi.vesicle.s,50);
%     xlabel('Vesicle amplitude s');
%     drawnow;
%
%     mi.basePath=ParsePath(inPath);  % make it the local path
%     bname=[mi.basePath mi.procPath mi.baseFilename];
%     WriteMRC(mv,ds0*mi.pixA,[bname 'mv.mrc']);
%     jname=[mi.basePath mi.procPath 'jpeg/' mi.baseFilename];
%     imwrite(uint8(imscale(rot90(mv),256,1e-3)),[jname 'mv.jpg']);
%
