function [RATES_corr] = gyro_correction(RATES, n, S, T, M)
% Correction of the measured angular rates using estimated sensor error model
% Inputs:
% RATES - measured angular rates
% n - number of measured data samples
% S - scale factor matrix
% T - orthogonalization matrix
% M - orthogonalization matrix

% Output:
% RATES_corr - corrected angular rates using sensor error model

for i = 1:n
        RATES_corr(i,:) = M^-1*T^-1*S^-1*RATES(i,:)';           % Correction of angular rates using estimated SEM
end