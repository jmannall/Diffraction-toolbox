close all

filepath = ['C:' filesep 'Documents' filesep 'GitHub' filesep 'DiffractionPlugin' filesep 'UnitTestData'];
CheckFileDir(filepath);

controlparameters = struct('fs', 48000, 'nfft', 4096, 'noDirect', true, 'saveFiles', 0);

zW = [2.5 10 7 5 2];
tW = [270 320 350 212 290];
tS = [45 10 2 60 33];
tR = [226 318 207 198 200];
rS = [1 2 7 3.5 2.7];
rR = [1 3.1 4 2 1.8];
zS = [1 3.4 3.2 2 1.7];
zR = [1 4.2 2.1 2.5 1.9];

out = [zW; tW; tS; tR; rS; rR; zS; zR];

filename = 'diffractionPaths.csv';
writematrix(out, [filepath filesep filename]);

out = [];
for i = 1:length(zW)
    ir = SingleWedge(zW(i), tW(i), tS(i), tR(i), rS(i), rR(i), zS(i), zR(i), controlparameters, true);
    idx = find(ir.diff1, 1, 'first');
    out = [out; ir.diff1(idx:end - 1)];
end

filename = 'btm.csv';
filepath = ['C:' filesep 'Documents' filesep 'GitHub' filesep 'DiffractionPlugin' filesep 'UnitTestData'];
CheckFileDir(filepath);

writematrix(out, [filepath filesep filename]);