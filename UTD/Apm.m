function [output, bArg] = Apm(n, B, plus)
    
    if plus
        N = round((pi + B) / (2 * pi * n));
    else
        N = round((-pi + B) / (2 * pi * n));
    end
    bArg = B - 2 * n * pi * N;
    output = 2 * cos(bArg / 2).^2;
    bArg = PlusMinus(bArg, plus);
end