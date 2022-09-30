close all
controlparameters = struct('fs',48000, 'difforder', 1, 'nfft', 8192);
controlparameters.difforder = 3;

x = SinglePanel(0.5, 2, 2, controlparameters, true);