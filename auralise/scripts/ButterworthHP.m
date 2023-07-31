function H = ButterworthHP(f, fc)
    s = 1j * f / fc;
    H = s .^ 2 ./ (1 + sqrt(2) * s + s .^ 2);
end