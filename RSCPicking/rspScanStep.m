function [scan,dis,done]=rspScanStep(mode,scan,dis,nFound)
% function [scan,dis,done]=rspScanStep('next',scan,dis,nFound)
%                     scan=rspScanStep('init',scan,dis); % initialization
                    
% handle looping for parameter scan in SimpleRSPicker.

if  strcmpi(mode,'init')
    scan.active=(dis.roboFitStep>0);
    scan.iAmp=double(scan.active);
    scan.iVar=0;
    scan.nAmps=inf;
    scan.ampStep=1.02;
    scan.varStep=.95;
    scan.amps=zeros(0,1);
    scan.vars=zeros(0,1);
    scan.nFound=[];
    if scan.active
        scan.basePars=dis.pars;
    elseif ~scan.active
        dis.pars=scan.basePars;
    end;
    return
elseif ~strcmpi(mode,'next')
    error(['Unrecognized mode: ' mode]);
end;

done=false;
if ~scan.active % we're not scanning
    return
end;
%   --------scanning--------
if scan.iVar==0 % first time through, no data yet, use the original pars
    scan.iVar=1;
    scan.nAmps=inf;
    return
else
    % update the statistics
    scan.amps(scan.iAmp)=dis.pars(1);
    scan.nParticles(scan.iAmp,scan.iVar)=nFound;
    if scan.iVar==1
        scan.amps(scan.iAmp)=dis.pars(1);
        if nFound==0 % Gotten to highest threshold
            scan.nAmps=scan.iAmp; % last iAmp value
            disp([num2str(scan.nAmps) ' amplitude vales being scanned.']);
        end;
    end;
%     if we're done, save the data
    if scan.iAmp==1 % check at the beginning of a scan
        scan.vars(scan.iVar)=dis.pars(10);
        if nFound==0 % we're all done, var is so low that even at the
%             first amp threshold we find no particles.
            scan.infoPath=dis.infoPath;
            scan.infoName=dis.infoName;
            baseName=scan.infoName(1:end-6);
            CheckAndMakeDir('Picker/',1);
            save(['Picker/ ' baseName 'ps.mat'],'scan');
            done=true;
        end;
    end;
    
    if ~done % update the parameters for the next search
        scan.iAmp=scan.iAmp+1;
        if scan.iAmp>scan.nAmps
            scan.iAmp=1;
            scan.iVar=scan.iVar+1;
        end;
        dis.pars(1)=scan.basePars(1)*scan.ampStep^(scan.iAmp-1);
        dis.pars(10)=scan.basePars(10)*scan.varStep^(scan.iVar-1);
    end;
end;
