function [output, bArg] = Apm(n, B, plus) % 13
    
    if plus
        N = round((pi + B) / (2 * pi * n)); % 5
    else
        N = round((-pi + B) / (2 * pi * n));
    end
    bArg = 2 * n * pi * N - B; % 4
    output = 2 * cos(bArg / 2).^2; % 4
    bArg = PlusMinus(-bArg, plus);
end