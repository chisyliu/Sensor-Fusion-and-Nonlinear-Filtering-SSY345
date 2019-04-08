clearvars; close all; clc;


%% Test genLinearStateSequence

% Tolerance
tol = 1e-1;

N = 5000;

% Define prior
x_0     = [0;0]; 
n       = length(x_0); 
P_0     = diag(ones(n,1));

% Define process model
A       = [1 1; 0 1];
Q       = diag(ones(n,1));

% generate state sequence
s = rng;
X = genLinearStateSequence(x_0, P_0, A, Q, N);

% Plot results
figure(1);clf;hold on;
plot(X');
title('Your solution');
xlabel('k');
ylabel('x');


%% Test genLinearMeasurementSequence
close all;

absTol = 1e-1;
relTol = 5e-2;

N = 50000;

n = randi(1,1);
m = randi(n,1);

% Define state sequence
X = zeros(n,N+1);

% Define measurement model
H = 1;
R = .5^2;

% Generate measurements
Y = genLinearMeasurementSequence(X, H, R);

assert(size(Y,1) == m, 'Y has the wrong measurement dimension');
assert(size(Y,2) == N, 'Y should have N columns');

% PLot results
figure(1);clf;hold on;
plot(0:10,X(1,1:11), '--k');
plot(1:10, Y(1,1:10), '*r');
legend('State sequence', 'Measurements')
title('Your solution');
xlabel('k');
ylabel('position');

%% Test kalmanFilter



% General parameters
N = 10;
n = 1;
T = 1;

% Motion Parameters
A = 1;
Q = T*.5^2;

% Measurement Parameters
H = 1;
R = 1^2;

% Prior
x_0  = 0;
P_0  = 2^2;

% genereate measurement sequence
Y = 10*ones(1,N);

% Filter
[stateSequence, covarianceSequence] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

