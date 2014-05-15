function [ANGLES] = gyro_integration(RATES, INTEGRATION_start, INTEGRATION_end, Ts)
% Integration of angular rate to angle for all axes
% Note: this integration method assumes rotation only around the sensitive
% axis of the gyroscope being calibrated
%
% RATES - Measured angular rates for integration
% INTEGRATION_start - first sample of angular rates for integration
% INTEGRATION_end - last sample of angular rate being integrated
% Ts - sampling time

ANGLES(1:3) = 0;
for i = INTEGRATION_start:INTEGRATION_end
    ANGLES = ANGLES + RATES(i,:)*Ts;
end
return