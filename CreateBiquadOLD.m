function [tfmag, fvec] = CreateBiquadOLD(zR, zI, pR, pI, k)
    
    fs = 48000;
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
    x = fft(b, 8192);
    y = fft(a, 8192);
    
    x = squeeze(prod(x,2));
    y = squeeze(prod(y,2));
    
    tfcomplex = x ./ y;
    epsilon = 1e-8;
    tfmag = 20*log(abs(tfcomplex(1:(end / 2),:))+epsilon) / log(10);
    




    z = [zR + zI * 1i; zR - zI * 1i];
    p = [pR + pI * 1i; pR - pI * 1i];

    numFreq = 4096;
    [bref, aref] = deal(zeros(numObservations, 2 * numBiquads + 1));
    href = zeros(numFreq, numObservations);
    z = extractdata(z);
    p = extractdata(p);
    k = extractdata(k);
    for i = 1:numObservations
        [bref(i,:), aref(i,:)] = zp2tf(z(:,i), p(:,i), k(i));
        [href(:,i), f] = freqz(bref(i,:), aref(i,:), numFreq);
    end
    
    fvec = (f * fs) / (2 * pi);
    tfmagref = 20*log10(abs(href));

%     f3 = figure(3);
%     clf(f3)
%     semilogx(fvec, tfmagref)
%     hold on
%     semilogx(fvec, extractdata(tfmag), '--')
%     hold off
end

