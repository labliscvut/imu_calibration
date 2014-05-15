function [ACC_corr] = acc_correction(ACC, n, T, S, b)
% Correction of the measured data using estimated sensor error model
% Inputs:
% data - measured accelerations
% n - number of measured orientation
% T - orthogonalization matrix
% S - scale factor matrix
% b - vector of offsets
% Output:
% ACC_corr - corrected accelerations using sensor error model

for i = 1:n
    ACC_corr(:,i) = T*S*(ACC(:,i) - b);             % Correction of accelerations using SEM estimated using fminunc function       
end
