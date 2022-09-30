close all
clear all

fs = 48000;

% CreateIIR(-0.5, 0.5, 0.5, @UpdateFirstOrderIIR, fs)
% CreateIIR(0.8, 0.4, 0.05, @UpdateFirstOrderIIR, fs)
% 
% CreateIIR([-0.9; 0.8], [0.85; 0.4], 0.05, @UpdateSecondOrderIIR, fs)

S.z = 0.8;
S.p = 0.4;
S.k = 0.05;

plotFcn = @(S) IIRFilter(S);
sliderFcn = @(S) CreateFirstOrderIIRSliders(S);
updateFcn = @UpdateFirstOrderIIR;

CreateInteractivePlot(S, plotFcn, sliderFcn, updateFcn, fs);

S.z = [-0.9; 0.8];
S.p = [0.85; 0.4];
S.k = 0.05;

plotFcn = @(S) IIRFilter(S);
sliderFcn = @(S) CreateSecondOrderIIRSliders(S);
updateFcn = @UpdateSecondOrderIIR;

CreateInteractivePlot(S, plotFcn, sliderFcn, updateFcn, fs);

%% Notes
% Hsh (zero > pole) gain between DC and fs / 2 determined by distance between zero and
% pole. Location of fc determied by average position between DC and fs / 2.

% Lpf zero at fs / 2. (As zero gets increased becomes a Lsh). Gain increased by
% moving pole towards DC. fc decreased as pole moved towards DC.

% Lpf: DC gain = 20log10(S.k * ((1-S.z) / (1-S.p)))     
%      fc = (y * (S.p^2 + 1) - S.z^2 - 1) / (2 * (S.p * y - S.z))
%       where y = ((1-S.z) ^ 2) / ((1-S.p) ^ 2 * 10 ^ 0.3)
% When S.z == -1
%      DC gain = 20log10(2 * S.k / (1 - S.p))
%      fc = (y * (S.p^2 + 1) - 2) / (2 * (S.p * y + 1));
%       where y = 4 / ((1 - S.p) ^ 2 * 10 ^ 0.3)

% Hsh: nyq gain = 20log10(S.k * ((1+S.z) / (1+S.p)))
%      B0 = DC gain / nyq gain
%      B0 = 20log10((1 - S.p) * (1 + S.z) / ((1 + S.p) * (1 - S.z)))

% For a give BTM with -5 20 gain, -20 20000 20000
% B0 = 











