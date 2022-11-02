function y = ConvolveTfcomplex(x, tfcomplex, n)
    H = ([tfcomplex; flipud(conj(tfcomplex(2:end,:)))])';
    h = ifft(H);
    h = h(1:round(length(H) / 2));
    % Code to perform Convolution using Overlap Add Method
    n1 = length(x);
    n2 = length(h);
    N = n1+n2-1;
    y = zeros(1,N);
    h1 = [h zeros(1,n2-1)];
    n3 = length(h1);
    y = zeros(1,N+n3-n2);
    
    numPerPos = floor((n1 - n2) / n);
    for i = 1:n2:n1 - n2
        if i<=(n1-n2-1)
            x1 = [x(i:i+n3-n2) zeros(1,n3-n2)];
        else
            x1 = [x(i:n1) zeros(1,n3-n2)];
        end
        x2 = fft(x1);
        x3 = x2.*H(ceil(i / numPerPos),:);
        x4 = ifft(x3);
        if (i==1)
            y(1:n3) = x4(1:n3);
        else
            y(i:i+n3-1) = y(i:i+n3-1)+x4(1:n3);
        end
    end
end