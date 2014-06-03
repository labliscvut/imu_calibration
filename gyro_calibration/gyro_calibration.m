%% SCRIPT FOR 3-AXIS GYROSCOPE CALIBRATION
% Calibration procedure: 
%   1) DATA PREPROCESSING: offset interval is determined for gyroscope offset estimation; intervals of the angular rates are marked for
%   all 3-rotations around the sensitive axes (x,y,z) (INTEGRATION_x_start, INTEGRATION_x_end, etc.)
%   2) INTEGRATION OF ANGULAR RATES
%   3) CALIBRATION ALGORITHM
%   4) CORRECTION OF THE MEASURED DATA: _rotation are corrected by the estimated SEM
%   5) EVALUATION
% 
% Authors: M. Sipos, J. Rohac, J. Simanek, Department of Measurement, Faculty of Electrical Engineering, Czech Technical University in Prague, Czech Republic
%   [1] Rohac, J.; Sipos, M.; Simanek, J.: Calibration of the Low-cost Triaxial Inertial Sensors


clc; clear all; close all; warning off;
fprintf('Triaxial gyroscope calibration \n\n');

load('measured_angular_rates.mat')                                  % Measured angular rates, rotations around x, y, z, axes

Ts = 0.01;                                                          % Sampling time
n = length(x_rotation);                                             % Length of measured data

ANGLES_ref = [361.338, 367.160, 363.6703];                          % Referential angles of rotations

OFFSET_start = 1;                                                   % Data interval - Angular rates under static conditions
OFFSET_end = 2000;                                                  % for estimation of offsets
 
INTEGRATION_x_start = 2420;                                         % Interval for integration of angular rates to obtain angle of rotation in x axis
INTEGRATION_x_end = 2910;                                           % Interval is determined as angular rate higher than threshold based on static data
INTEGRATION_y_start = 2565;                                         % Interval for integration of angular rates to obtain angle of rotation in y axis
INTEGRATION_y_end = 3120;                                           % Interval is determined as angular rate higher than threshold based on static data
INTEGRATION_z_start = 2975;                                         % Interval for integration of angular rates to obtain angle of rotation in z axis
INTEGRATION_z_end = 3500;                                           % Interval is determined as angular rate higher than threshold based on static data
    
OFFSET_rot_x = mean(x_rotation(OFFSET_start:OFFSET_end,:));         % Computation of offset from rotation around x axis
OFFSET_rot_y = mean(y_rotation(OFFSET_start:OFFSET_end,:));         % Computation of offset from rotation around y axis
OFFSET_rot_z = mean(z_rotation(OFFSET_start:OFFSET_end,:));         % Computation of offset from rotation around z axis
OFFSET = [OFFSET_rot_x(1), OFFSET_rot_x(2), OFFSET_rot_x(3)];       % Vector of offsets

x_rotation = x_rotation - repmat(OFFSET_rot_x, n, 1);               % Substraction of offset from angular rates
y_rotation = y_rotation - repmat(OFFSET_rot_y, n, 1);               % Substraction of offset from angular rates
z_rotation = z_rotation - repmat(OFFSET_rot_z, n, 1);               % Substraction of offset from angular rates

%% --- Interation of angular rates ---
[Yg(1,:)] = gyro_integration(x_rotation, INTEGRATION_x_start, INTEGRATION_x_end, Ts);   % Computation of angles of rotation - (3x3) Yg matrix - x axis
[Yg(2,:)] = gyro_integration(y_rotation, INTEGRATION_y_start, INTEGRATION_y_end, Ts);   % Computation of angles of rotation - (3x3) Yg matrix - y axis
[Yg(3,:)] = gyro_integration(z_rotation, INTEGRATION_z_start, INTEGRATION_z_end, Ts);   % Computation of angles of rotation - (3x3) Yg matrix - z axis

%% --- Calibration algorithm ---
Ug = diag([ANGLES_ref(1),ANGLES_ref(2),ANGLES_ref(3)]);             % Diagonal reference angle matrix Ug
ST = chol((Yg*Ug^-1)*(Yg*Ug^-1)')';                                 % Cholesky decomposition -> multiplication of lower and upper triangular matrix - eq 8. [1]
[T,S] = lu(ST);                                                     % LU factorization -> lower and upper triangular matrix - eq. 9 [1]
M = T^-1*S^-1*Yg*Ug^-1;                                             % Computation of alignment matrix - eq. 10 [1]

%% --- Correction of measured data ---
[x_rotation_corr] = gyro_correction(x_rotation, n, S, T, M);        % Correction of measured angular rates using estimated SEM 
[y_rotation_corr] = gyro_correction(y_rotation, n, S, T, M);        % (measured angular rates, number of samples, T, SF, b)
[z_rotation_corr] = gyro_correction(z_rotation, n, S, T, M);

% constructs the Yg matrix for corrected angular rates
[Yg_corr(1,:)] = gyro_integration(x_rotation_corr, INTEGRATION_x_start, INTEGRATION_x_end, Ts); % Computation of angles of rotation - (3x3) Yg matrix - x axis
[Yg_corr(2,:)] = gyro_integration(y_rotation_corr, INTEGRATION_y_start, INTEGRATION_y_end, Ts); % Computation of angles of rotation - (3x3) Yg matrix - y axis
[Yg_corr(3,:)] = gyro_integration(z_rotation_corr, INTEGRATION_z_start, INTEGRATION_z_end, Ts); % Computation of angles of rotation - (3x3) Yg matrix - z axis

%% --- Results of angular rates calibration ---
fprintf('SENSOR ERROR MODEL:\n\n');
fprintf('Vector of offsets (deg/s):\n'); disp(OFFSET);
fprintf('Scale-factor matrix (-):\n'); disp(S);
fprintf('Orthogonalization matrix (-):\n'); disp(T);
fprintf('Alignment matrix (-):\n'); disp(M);

fprintf('EVALUATION OF CORRECTION:\n\n');
fprintf('Yg matrix (deg):\n'); 
fprintf(' rotation x projected to x, rotation x projected to y, rotation x projected to z\n rotation y projected to x, rotation y projected to y, rotation y projected to z\n rotation z projected to x, rotation z projected to y, rotation z projected to z\n');
fprintf('Reference angles of rotation (deg):\n'); disp(Ug);
fprintf('Angles of rotation obtained by integrating angular rates before correction (deg):\n'); disp(Yg);
fprintf('Angles of rotation obtained by integrating angular rates after correction (deg):\n'); disp(Yg_corr);
fprintf('Deviations from reference angles - Angles of rotation obtained by integrating angular rates before calibration (deg):\n'); disp(Yg-Ug);
fprintf('Deviations from reference angles - Angles of rotation corrected by deterministic sensor error model (deg): \n'); disp(Yg_corr-Ug);

figure; 
plot((1:n)*Ts, x_rotation); grid on;
title('Measured angular rates - x axis');
xlabel('sampling time (s)'); ylabel('measured angular rate (deg/s)'); legend('x axis', 'y axis', 'z axis');
figure; 
plot((1:n)*Ts, y_rotation);  grid on;
title('Measured angular rates - y axis');
xlabel('sampling time (s)'); ylabel('measured angular rate (deg/s)'); legend('x axis', 'y axis', 'z axis');
figure; 
plot((1:n)*Ts, z_rotation);  grid on;
title('Measured angular rates - y axis');
xlabel('sampling time (s)'); ylabel('measured angular rate (deg/s)'); legend('x axis', 'y axis', 'z axis');

figure(111);
close(111);