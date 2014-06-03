%% SCRIPT FOR 3-AXIS ACCELEROMETER CALIBRATION
% Calibration procedure: 
%   1) DATA PREPROCESSING: 3-axis accelerometer output is measured and averaged in different static orientations (ACC_meas)
%   2) OPTIMIZATION: fminunc is used for Sensor Error Model (SEM) parameter estimation
%   3) EVALUATION: ACC_meas are corrected by the estimated SEM, root mean squared error (RMSE) is computed
% 
% Authors: M. Sipos, J. Rohac, J. Simanek, Department of Measurement, Faculty of Electrical Engineering, Czech Technical University in Prague, Czech Republic
%   [1] Rohac, J.; Sipos, M.; Simanek, J.: Calibration of the Low-cost Triaxial Inertial Sensors

clc; close all; clear all; warning off;
fprintf('Triaxial accelerometer calibration based on FMINUNC function \n');

load('measured_accelerations.mat');        % Data for calibration - each vector ACC_meas(1:3,i) is a mean of 1000 static samples @ 100 Hz (g)                     
G = 1;                                     % Magnitude of the gravity field vector (~9.81 m/s^2)
n = length(ACC_meas);                      % Number of samples in ACC_meas

x0 = [0, 0, 0, 1, 1, 1, 0, 0, 0];          % Initial values of nine-state vector using FMINUNC
                                           % Structure of vector: x0 = [ALFAxy, ALFAzx, ALFAzy, Sx, Sy, Sz, bx, by, bz];

[x, fval, exitflag, output, grad] = fminunc(@(x)acc_minimization(x, n, ACC_meas, G), x0);   % Minimization using fminunc function

fprintf('Sensor Error Model: \n');
% Sensor Error Model 
T = [1, 0, 0; x(1), 1, 0; x(2), x(3), 1]                % Orthogonalization matrix
S = [x(4), 0, 0; 0, x(5), 0; 0, 0, x(6)]                % Scale Factor Matrix
b = [x(7); x(8); x(9)]                                  % Offsets

[ACC_corr] = acc_correction(ACC_meas, n, T, S, b);      % Correction of accelerations using SEM estimated using fminunc function 
                                                        % (measured_data, number of samples, T, S, b)

[RMSE(1), delta(1,:)] = acc_rmse(ACC_meas, n, G);       % Computation of RMSE before calibration
[RMSE(2), delta(2,:)] = acc_rmse(ACC_corr, n, G);       % Computation of RMSE after calibration

fprintf('RMSE before calibration:     %.6f g \n', RMSE(1));
fprintf('RMSE after calibration:      %.6f g \n', RMSE(2));

figure; hold on; grid on;
plot(delta(1,:), 'b--');
plot(delta(2,:), 'r');
title('Deviations \Delta before and after correction using SEM in measured orienations'); 
xlabel('measured orientation (-)'); 
ylabel('\Delta = sqrt(acc_x^2+acc_x^2+acc_z^2) - G (g)');
legend('before correction','after correction');