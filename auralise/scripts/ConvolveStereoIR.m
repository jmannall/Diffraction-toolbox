function audio = ConvolveStereoIR(input, ir, windowLength)
    audio.L = ConvolveIR(input, ir.L, windowLength);
    audio.R = ConvolveIR(input, ir.R, windowLength);
end