%% Calculate the transfer function and magnitude response of an impulse response

function [tfmag, tfcomplex, ir] = IrToTf(ir, nfft)  

    structInput = isstruct(ir);

    if structInput
        fields = fieldnames(ir);
        numFields = length(fields);
        numReceivers = size(ir.(fields{1}), 2);
        for n = 1:numFields
            field = fields{n};

            F = fft(ir.(field),nfft);
            tfcomplex.(field) = F(1:nfft/2,:);
            tfmag.(field) = 20*log10(max(abs(F(1:nfft/2,:)), 1e-10));
        end
    else
        F = fft(ir,nfft);
        tfcomplex = F(1:nfft/2,:);
        tfmag = mag2db((max(abs(tfcomplex), 1e-10)));
    end
end