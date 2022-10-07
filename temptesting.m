
fs = 48000;

pathlength = 1;
c = 344;
delay = pathlength / c;
sampledelay = delay * fs;

pathlength = round(sampledelay) / fs * c;

frac = 1.2;
step = c * (1 / fs);

delay = 2 * pathlength / c;
sampledelay1 = delay * fs + 1;
delay = (2 * pathlength + step) / c;
sampledelay2 = delay * fs + 1;
delay = (2 * pathlength + frac * step) / c;
sampledelay3 = delay * fs + 1;

controlparameters = struct('difforder', 1, 'nfft', 8192', 'fs', fs);

[ir1, tfmag1] = SingleWedge(20, 350, 10, 320, pathlength, pathlength, 10, 10, controlparameters, false);
[ir2, tfmag2] = SingleWedge(20, 350, 10, 320, pathlength, pathlength + step, 10, 10, controlparameters, false);
[ir3, tfmag3] = SingleWedge(20, 350, 10, 320, pathlength, pathlength + frac * step, 10, 10, controlparameters, false);

ir1 = ir1.diff1;
ir2 = ir2.diff1;
ir3 = ir3.diff1;