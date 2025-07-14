function [b, a] = PeakFilters(fc, g, Q, fs)

    %fc = fc(2:end - 1);
    %g = g(2:end - 1);
    
    omega = 2 * pi * fc(2:end - 1) / fs;
    alpha = sin(omega) / (2 * Q);

    r = 1 ./ (1 + alpha);

    b = r .* [1 + g(2:end -1) .* alpha; -2 * cos(omega); 1 - g(2:end -1) .* alpha];
    a = [ones(1,length(fc) - 2); b(2,:); r .* (1 - alpha)];


    alpha = sin(omega) / ( 2 * Q);
    A = sqrt(g(2:end - 1)); % 10 .^ (gdb / 40);

    v1 = alpha .* A;
    v2 = alpha ./ A;
    v3 = -2 * cos(omega);

    b2 = [1 + v1; v3; 1 - v1];
    a2 = [1 + v2; v3; 1 - v2];

    b = b2 ./ a2(1,:);
    a = a2 ./ a2(1,:);

    omega = 2 * pi * fc(1) / fs;
    alpha = sin(omega) / ( 2 * Q);
    A = sqrt(g(1)); % 10 .^ (gdb / 40);

    v1 = (A + 1);
    v2 = (A - 1);
    v3 = v1 .* cos(omega);
    v4 = v2 .* cos(omega);
    v5 = 2 * sqrt(A) .* alpha;

    alphaTest = sin(omega) / Q;
    v5Test = sqrt(A) .* alphaTest;

    bLowShelf = A .* [ v1 - v4 + v5; 2 * (v2 - v3); v1 - v4 - v5];
    aLowShelf = [v1 + v4 + v5; -2 * (v2 + v3); v1 + v4 - v5];

    bLowShelf = bLowShelf ./ aLowShelf(1,:);
    aLowShelf = aLowShelf ./ aLowShelf(1,:);

    omega = 2 * pi * fc(end) / fs;
    alpha = sin(omega) / ( 2 * Q);
    A = sqrt(g(end)); % 10 .^ (gdb / 40);

    v1 = (A + 1);
    v2 = (A - 1);
    v3 = v1 .* cos(omega);
    v4 = v2 .* cos(omega);
    v5 = 2 * sqrt(A) .* alpha;

    alphaTest = sin(omega) / Q;
    v5Test = sqrt(A) .* alphaTest;

    bHighShelf = A .* [ v1 + v4 + v5; -2 * (v2 + v3); v1 + v4 - v5];
    aHighShelf = [v1 - v4 + v5; 2 * (v2 - v3); v1 - v4 - v5];

    bHighShelf = bHighShelf ./ aHighShelf(1,:);
    aHighShelf = aHighShelf ./ aHighShelf(1,:);

    b = [bLowShelf, b, bHighShelf];
    a = [aLowShelf, a, aHighShelf];
end