function [RMSE, delta] = acc_rmse(ACC, n, G)
% Computes the Root Mean Squared Error (RMSE) and the deviation (dev) between the norm of ACC and magnitude of gravity acceleration
% Inputs:
% n - number of measured orientation
% data - measured accelerations
% G - magnitude of gravity acceleration

% Outputs:
% RMSE - Root Mean Squared Error
% dev  - deviation of sqrt(acc_x^2+acc_x^2+acc_z^2)-G

dRMSE = 0;
delta = zeros(1,n);
for i = 1:n
    delta(1,i) = sqrt(ACC(1,i)^2 + ACC(2,i)^2 + ACC(3,i)^2) - G;   
    dRMSE = dRMSE + delta(i)^2;       
end
RMSE = sqrt(dRMSE/n);
