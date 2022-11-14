function y = FftConvolution(x, H)
    Nx = length(x);
    Nh = length(ifft(H));
    Ny = Nx + Nh - 1;

    K = nextpow2(Ny);
    X = fft([x; zeros(2 ^ K - Nx, 1)]);
    H = fft([ifft(H); zeros(2 ^ K - Nh, 1)]);

    Y = X .* H;
    y = ifft(Y);
end