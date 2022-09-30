clear all


zR = dlarray([0.5; 0.1; -0.2]);
zI = dlarray([0.2; 0.3; 0.1]);
pR = dlarray([-0.2; -0.8; 0.1]);
pI = dlarray([0.1; 0.4; 0.2]);
k = dlarray(0.3);

tfvalue = zeros(4096, 1);

for i = 1:8000
    [fval, gradval] = dlfeval(@GradientTest, zR, zI, pR, pI, k, tfvalue);

    zR = zR - gradval{1} * 0.1;
    zI = zI - gradval{2} * 0.1;
    pR = pR - gradval{3} * 0.1;
    pI = pI - gradval{4} * 0.1;
    k = k - gradval{5} * 0.1;

    tfmag = CreateBiquad(zR, zI, pR, pI, k);
    if mod(i, 500) == 0
        disp(['Iteration: ', num2str(i)]);
        disp(['Loss: ', num2str(fval)]);
    end

end

%% Data
z = [zR + zI * 1i; zR - zI * 1i];
p = [pR + pI * 1i; pR - pI * 1i];

z = extractdata(z);
p = extractdata(p);
k = extractdata(k);

epsilon = 1e-8;

[bold, aold] = zp2tf(z, p, k);

z = (1 - epsilon) .* z .* tanh(abs(z)) ./ abs(z + epsilon);
p = (1 - epsilon) .* p .* tanh(abs(p)) ./ abs(p + epsilon);

[b, a] = zp2tf(z, p, k);

[hold_, ~] = freqz(bold, aold, 2048);
[h, f] = freqz(b, a, 2048);

figure
semilogx(f, 20*log10(abs(h)))
hold on
semilogx(f, 20*log10(abs(hold_)))



function [y, dydx] = GradientTest(zR, zI, pR, pI, k, tfvalue)

    tfmag = CreateBiquad(zR, zI, pR, pI, k);
    y = sum((tfmag - tfvalue).^2) / 4096;
    dydx  = dlgradient(y, {zR, zI, pR, pI, k});

end