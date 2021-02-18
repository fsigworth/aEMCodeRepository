% rlDisplayModelSpectra

starName='run_model.star';

[nm,da]=ReadStarFile(starName);
gpString='data_model_group_';
nGpString=numel(gpString);

nnm=numel(nm);

for ind=1:nnm
    if strcmp(nm{ind},'data_model_groups')
        dGroups=da{ind};
    elseif strncmp(nm{ind},gpString,nGpString);
        d=da{ind};
        str=nm{ind}(nGpString+1:end);
        groupInd=str2double(str);
%          groupInd=str2double(nm{ind}(nGpString+1,end)); % Matlab bug?
%          groupInd=str2num(nm{ind}(nGpString+1,end)); % Matlab bug?
        if groupInd==1
            spectra=d.rlnSigma2Noise;
            freqs=d.rlnResolution;
        else
            spectra(:,groupInd)=d.rlnSigma2Noise;
        end;
    end;
end;
numSpectra=size(spectra,2)

subplot(211);
plot(dGroups.rlnGroupScaleCorrection);
xlabel('Group');
ylabel('Scale Correction');

subplot(212);
semilogy(freqs,spectra);
xlabel('Frequency, A^-1');
ylabel('Spectral density');

