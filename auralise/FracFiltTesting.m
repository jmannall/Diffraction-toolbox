close all
clear all
set(0, 'DefaultLineLineWidth', 2);

%% Create frac filters

fs = 48000;

len = [4 8 16];
fd = 0.2;

[h1, i10, bw1] = designFracDelayFIR(fd,len(1));
[h2, i20, bw2] = designFracDelayFIR(fd,len(2));
[h3, i30, bw3] = designFracDelayFIR(fd,len(3));

x1 = (-len(1):2:len(1) - 1) / 2;
x2 = (-len(2):2:len(2) - 1) / 2;
x3 = (-len(3):2:len(3) - 1) / 2;
p11 = hamming(len(1))' .* sinc(x1 - fd + 1);
p21 = hamming(len(2))' .* sinc(x2 - fd + 1);
p31 = hamming(len(3))' .* sinc(x3 - fd + 1);
p12 = hanning(len(1))' .* sinc(x1 - fd + 1);
p22 = hanning(len(2))' .* sinc(x2 - fd + 1);
p32 = hanning(len(3))' .* sinc(x3 - fd + 1);

%% Stems

figure
stem(x1, h1, '-o')
hold on
stem(x1, p11, '-o')
stem(x1, p12, '-o')
legend('MATLAB', 'Sinc hamming', 'Sinc hanning')

figure
stem(x2, h2, '-o')
hold on
stem(x2, p21, '-o')
stem(x2, p22, '-o')
legend('MATLAB', 'Sinc hamming', 'Sinc hanning')

figure
stem(x3, h3, '-o')
hold on
stem(x3, p31, '-o')
stem(x3, p32, '-o')
legend('MATLAB', 'Sinc hamming', 'Sinc hanning')

%% Magnitude

[H1, w] = freqz(h1, 1, [], fs);
H2 = freqz(h2, 1);
H3 = freqz(h3, 1);

P11 = freqz(p11, 1);
P21 = freqz(p21, 1);
P31 = freqz(p31, 1);
P12 = freqz(p12, 1);
P22 = freqz(p22, 1);
P32 = freqz(p32, 1);

figure
semilogx(w,mag2db(abs([H1 P11 P12])))
hold on
semilogx(w,mag2db(abs([H2 P21 P22])), '--')
semilogx(w,mag2db(abs([H3 P31 P32])), '-.')
legend('MATLAB', 'Hamming', 'Hanning', 'MATLAB', 'Hamming', 'Hanning', 'MATLAB', 'Hamming', 'Hanning', 'Location', 'southwest')

% Look at implementing farrow based larange fractional filters


fs = 48000;

delay = 0.5;
fd = rem(delay, 1);
len = [4 8 16];
[h1,i10,bw1] = designFracDelayFIR(fd,len(1));
[h2,i20,bw2] = designFracDelayFIR(fd,len(2));
[h3,i30,bw3] = designFracDelayFIR(fd,len(3));

[H1,w] = freqz(h1,1,[],fs);
H2 = freqz(h2,1);
H3 = freqz(h3,1);

v1 = ones(1, len(1));
v2 = ones(1, len(2));
v3 = ones(1, len(3));
for i = 2:len(1)
    v1(i) = ((2 * delay + 1) / (i - 1) - 1) * v1(i-1);
end
for i = 2:len(2)
    v2(i) = ((2 * delay + 1) / (i - 1) - 1) * v2(i-1);
end
for i = 2:len(3)
    v3(i) = ((2 * delay + 1) / (i - 1) - 1) * v3(i-1);
end

p1 = 4^(-delay) * v1;
p2 = 4^(-delay) * v2;
p3 = 4^(-delay) * v3;

P1 = freqz(p1,1);
P2 = freqz(p2,1);
P3 = freqz(p3,1);

%w = w * fs / (2 * pi);

figure;
semilogx(w,mag2db(abs([H1 H2 H3])))
hold on
semilogx(w,mag2db(abs([P1])))
semilogx(w,mag2db(abs([P2])), '--')
semilogx(w,mag2db(abs([P3])), '-.')
%hline = lines;
%xline(bw1, LineStyle = '--', LineWidth = 2, Color = hline(1,:))
%xline(bw2, LineStyle = '--', LineWidth = 2, Color = hline(2,:))
%xline(bw3, LineStyle = '--', LineWidth = 2, Color = hline(3,:))
hold off
title('Magnitude Responses in dB', FontSize = 12)
xlabel("Normalized Frequency (\times\pi rad/sec)")
ylabel("Magnitude (dB)")
grid
 
legend('Gain Response (N = 8)','Gain Response (N = 32)','Gain Response (N = 64)',...
       'Measured Combined Bandwidth (N = 8)',...
       'Measured Combined Bandwidth (N = 32)',...
       'Measured Combined Bandwidth (N = 64)', ...
       'Measured Combined Bandwidth (N = 100)',Location = 'Southwest')
ylim([-10 5])

[g1,w] = grpdelay(h1,1,[],fs);
g2 = grpdelay(h2,1);
g3 = grpdelay(h3,1);

d1 = grpdelay(p1,1);
d2 = grpdelay(p2,1);
d3 = grpdelay(p3,1);

%w = w * fs / (2 * pi);
figure;
semilogx(w,[g1 g2 g3])
hold on
semilogx(w,[d1 d2 d3])

hline = lines;
%yline(i10+fd, LineStyle = '--', LineWidth = 2, Color = hline(1,:))
%yline(i20+fd, LineStyle = '--', LineWidth = 2, Color = hline(2,:))
%yline(i30+fd, LineStyle = '--', LineWidth = 2, Color = hline(3,:))


 
title('Group Delay Responses',FontSize = 12)

xlabel("Normalized Frequency (\times\pi rad/sec)")
ylabel("Group Delay")
grid
 
legend('Group Delay Response (N = 8)','Group Delay Response (N = 32)','Group Delay Response (N = 64)',...
        'Nominal Group Delay (N = 8)',...
       'Nominal Group Delay (N = 32)',...
      'Nominal Group Delay (N = 64)', ...
      'Nominal Group Delay (N = 100)',Location = 'west',FontSize = 10);
ylim([-10,20]);

%% Function

fd = 0.3;
len = 28;

[h,i0,bw] = designFracDelayFIR(fd,len);

fdf = dsp.FIRFilter(h);

f = @(t) (0.1*t.^2+cos(0.9*t)).*exp(-0.1*(t-5).^2);
n = (0:38)'; t = linspace(0,38,512);
x = f(n); % Samples

scatter(n,x,20,'k','filled');
hold on;
plot(t,f(t),'color',[0.5 0.5 0.5],'LineWidth',0.5)
hold off;
xlabel('Time')
legend(["Input Samples","f(t)"])
title('Input Sequence with Known Underlying Analog Signal')
ax = gca; ax.XGrid='on';

y = fdf(x);

subplot(2,1,1);
stem(x);
title('Input Sequence');
xlabel('n')
subplot(2,1,2)
stem(y);
title('FIR Output Sequence');
xlabel('n')

figure
scatter(n,y,20,'red','filled')
hold on;
scatter(n+i0+fd,x,20,'black','filled')
plot(t,f(t-i0-fd),'Color',[1,0.5,0.5],'LineWidth',0.1)
xlabel('Time')
legend(["Filter output","Shifted Input Samples","Shifted f(t-i0-fd)"])
hold off
grid on
title('Input and Output Sequences Aligned and Overlaid')

%% Specified

fd  = 0.02;
tbw = 5/6;
[h,i0,mbw] = designFracDelayFIR(fd, tbw);

stem((0:length(h)-1),h); 
xlabel('h'); ylabel('h[n]');
title('Impulse Response of the Fractional Delay FIR')

[H1,w] = freqz(h,1);
G1 = grpdelay(h,1);
figure;
yyaxis left
plot(w/pi,mag2db(abs(H1)))
ylabel("Magnitude (dB)")

hold on
yyaxis right
plot(w/pi,G1)
ylabel("Group Delay (in samples)")

hline = lines;
xline(mbw, LineStyle =':',  Color = 'b', LineWidth = 2)
xline(tbw, LineStyle = '--', Color = 'm', LineWidth = 2)
yline(i0+fd, LineStyle = ':', Color = 'r', LineWidth = 1)
yticks([i0, i0+fd,i0+1:i0+9]); 

hold off
title('Magnitude Responses (dB) and Group Delay', FontSize = 10)
xlabel("Normalized Frequency (\times\pi rad/sample)")

legend('Gain Response','Group Delay Response','Measured Combined Bandwidth',...
        'Target Combined Bandwidth','Nominal Group Delay', ...
        Location = 'west', FontSize = 10)