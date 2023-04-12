function H = SingleUDFA(f, fc, g)

    alpha = 0.5;
    b = 1.44;
    Q = 0.2;
    r = 1.6;

    b = 1 + (b - 1) * g ^ 2;
    Q = 0.5 + (Q - 0.5) * g ^ 2;

    H = ((1j * f / fc) .^ (2 / b) + (1j * f / (Q * fc)) .^ (1 / (b ^ r)) + 1) .^ (-alpha * b / 2);
end