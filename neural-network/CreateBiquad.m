%% Create output of biquad filters from ZPK parameters. Traceable by dlgradient

function [tfmag, fvec, tfcomplex] = CreateBiquad(zR, zI, pR, pI, k, numFreq, fs)
    
    %% Calculate tfmag
    numObservations = size(zR,2);

    numBiquads = size(zR,1);
    
    [b, a] = deal(dlarray(zeros(3,numBiquads,numObservations)));
    for i = 1:numBiquads
        b(:,i,:) = [ones(1,numObservations); -2 * zR(i,:); zR(i,:).^2 + zI(i,:).^2];
        a(:,i,:) = [ones(1,numObservations); -2 * pR(i,:); pR(i,:).^2 + pI(i,:).^2];
    end
    b(:,1,:) = squeeze(b(:,1,:)) .* k;
    b = stripdims(b);
    a = stripdims(a);
    b = real(b);
    a = real(a);
    
    n = 2 * numFreq;
    x = fft(b, n);
    y = fft(a, n);
    
    x = squeeze(prod(x,2));
    y = squeeze(prod(y,2));
    
    F = x ./ y;
    tfcomplex = F(1:(end / 2),:);
    epsilon = 1e-8;
    tfmag = (20*log(abs(tfcomplex)+epsilon) / log(10))';

    %% Create fvec 

    if nargin == 6
        fvec = 0;
        disp('No sample rate provided for fvec');
    else
        [~, f] = freqz(b(:,1), a(:,1), numFreq);
        fvec = (f * fs) / (2 * pi);
    end
end

