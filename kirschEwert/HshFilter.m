function H = HshFilter(f, fc, g)

    s = (1j * f / fc)';
    H = (1 + sqrt(g) * s) ./ (1 + s / sqrt(g));
end