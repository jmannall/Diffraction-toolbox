function [output] = TanhFunctions(x, w, bA, mA)
    output = [tanh(1 ./ (x(1) * mA + x(2))), tanh(x(3) ./ w), tanh(x(4) ./ bA)];
end