function [zR, zI, pR, pI, k] = RandZPK(numBiquads, numFilters)
    Rz = sqrt(rand(numBiquads, numFilters));
    Argz = pi .* rand(numBiquads, numFilters);

    Rp = sqrt(rand(numBiquads, numFilters));
    Argp = pi .* rand(numBiquads, numFilters);

    zR = Rz .* cos(Argz);
    zI = Rz .* sin(Argz);
    pR = Rp .* cos(Argp);
    pI = Rp .* sin(Argp);

    k = 100 * rand(1, numFilters);
end