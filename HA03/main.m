close all; clear all; clc;

x = rand(5,1);
% Covariance
P = sqrt(10)*rand(5,5);
P = P*P';

% Random sensor position sequence
s1 = rand(2,1)*100;
s2 = rand(2,1)*100;

% Measurement model
h = @(x) dualBearingMeasurement(x, s1, s2);

% Random measurement
y = h(x)+rand(2,1);

% Measurement noise covariance
R = diag(rand(1,2).^2);

% EKF
[xp, Pp] = nonLinKFupdate(x, P, y, h, R, 'EKF')
