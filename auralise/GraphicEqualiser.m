function [tfmag, fvec, b, a] = GraphicEqualiser(material, fc)

fs = 48e3;
nfft = 4096;
Q = 0.98;

g = sqrt(1 - material);
gain = db2mag(mean(mag2db(g)));
target = mag2db(g);
g = g / gain;

fc = [(fc(1) / 2) * sqrt(fc(1) / (fc(1) / 2)), fc, fc(end) * sqrt((2.0 * fc(end)) / fc(end))];
g = [(g(1) + g(2)) / 2, g, (g(end - 1) + g(end)) / 2];

gdb = mag2db(g);
P = 6;

[b, a] = PeakFilters(fc, db2mag(P) * ones(size(fc)), Q, fs);

% [tfmag, fvec] = CalculateFilterResponse(b, a, nfft, fs);
% 
% o = 2 * pi * fvec / fs;
% num = (b(1,:) + b(2,:) .* exp(-1i * o) + b(3,:) .* exp(-1i * 2 * o));
% den = (a(1,:) + a(2,:) .* exp(-1i * o) + a(3,:) .* exp(-1i * 2 * o));
% 
% out = mag2db(abs(num ./ den));

%%
%f = [100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1e3, 1.25e3, 1.6e3, 2e3, 2.5e3, 3.15e3, 4e3, 5e3, 6.3e3, 8e3, 10e3];
%f = [200, 250, 315, 400, 500, 630, 800, 1e3, 1.25e3, 1.6e3, 2e3, 2.5e3, 3.15e3, 4e3, 5e3];
f = fc;
%f = [fc(1) / 2, fc, 2 * fc(end)];
f(1) = f(2) / 2;
f(end) = 2 * f(end - 1);
o = 2 * pi * f' / fs;
num = (b(1,:) + b(2,:) .* exp(-1i * o) + b(3,:) .* exp(-1i * 2 * o));
den = (a(1,:) + a(2,:) .* exp(-1i * o) + a(3,:) .* exp(-1i * 2 * o));

out = mag2db(abs(num ./ den));

if (length(f) ~= length(out))
    fm = zeros(1, length(fc) - 1);
    for i = 1:length(fc) - 1
        fm(i) = fc(i) * sqrt(fc(i + 1) / fc(i));
    end
    
    fidx = ones(size(f));
    
    for i = 1:length(fm)
        fidx(f > fm(i)) = i + 1;
    end
    
    for i = 1:length(fc)
        outBand(i,:) = mean(out(fidx == i, :), 1);
    end
    
    out = outBand;
end

out = out';
checkInv = inv(out);

inverse = eye(length(f));
for i = 1:length(f)
    factor = out(i,i);
    out(i,:) = out(i,:) / factor;
    inverse(i,:) = inverse(i,:) / factor;
    for j = 1:length(fc)
        if (j ~= i)
            factor = out(j,i) / out(i,i);
            out(j,:) = out(j,:) - factor * out(i,:);
            inverse(j,:) = inverse(j,:) - factor * inverse(i,:);
        end
    end
end

k = gdb * inverse;
gdb = k .* P;
g = db2mag(gdb);


%%
% o = 2 * pi * fc' / fs;
% num = (b(1,:) + b(2,:) .* exp(-1i * o) + b(3,:) .* exp(-1i * 2 * o));
% den = (a(1,:) + a(2,:) .* exp(-1i * o) + a(3,:) .* exp(-1i * 2 * o));
% 
% out = mag2db(abs(num ./ den));
% 
% k = gdb * inv(out');
% gdb = k .* P;
% g = db2mag(gdb);

%%

[b, a] = PeakFilters(fc, g, Q, fs);

[~, fvec] = CalculateFilterResponse(b, a, nfft, fs);

o = 2 * pi * fvec / fs;
num = (b(1,:) + b(2,:) .* exp(-1i * o) + b(3,:) .* exp(-1i * 2 * o));
den = (a(1,:) + a(2,:) .* exp(-1i * o) + a(3,:) .* exp(-1i * 2 * o));

tfmagI = mag2db(abs(num ./ den));

b(:,1) = gain * b(:,1);
[tfmag, fvec] = CalculateFilterResponse(b, a, nfft, fs);

semilogx(fvec, tfmag)
hold on
grid on
semilogx(fvec, tfmagI)
semilogx(fc(2:end - 1), target, 'x')
semilogx(fc, gdb, 'o')
end