function H = Butterworth(f, fc)
    s = 1j * f / fc;
    H = 1 ./ (1 + sqrt(2) * s + s .^ 2);
end