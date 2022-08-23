function output = Multiply(target, input)
    realPart = real(target) - real(input);
    imagPart = imag(target) - imag(input);
    distance = sqrt(realPart .^ 2 + imagPart .^ 2);
    output = prod(distance,"all");
end