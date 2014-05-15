function [RMSE] = acc_minimization(x, n, ACC, G)
% Minimization of criterion function - correction of data in each step and
% computation of RMSE as a criterion for minimization
%
% Inputs:
% x - vector of estimated parameters
% n - number of measured orientation
% data - measured accelerations
% G - magnitude of gravity acceleration
%
% Output:
% sigma - RMSE in actual step of minimization

T = [1,0,0;x(1),1,0;x(2),x(3),1];           % Sensor Error Model - orthogonalization matrix
S = [x(4),0,0;0,x(5),0;0,0,x(6)];           % SEM - Scale Factor Matrix
b = [x(7);x(8);x(9)];                       % SEM - Offsets

[ACC_est] = acc_correction(ACC, n, T, S, b);           % Correction of accelerations using SEM estimated using fminunc function 
[RMSE, delta] = acc_rmse(ACC_est, n, G);               % Computation of RMSE

return