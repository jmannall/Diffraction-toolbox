close all
clear all

fs = 48000;

CreateBTM(350, 1, 20, @UpdateBTM, fs)
% CreateBTM(320, 270, 10, @UpdateBTM, fs)
% 
% CreateIIR([-0.9; 0.8], [0.85; 0.4], 0.05, @UpdateSecondOrderIIR, fs)


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











