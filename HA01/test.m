clearvars; close all; clc;

%% Test sigmaElipse2D

%Here we give parameters for a Gaussian density. The parameter mu is the mean, and P is the covariance matrix.
mu = [-2; 1];
P = [4, -2; -2 2];

%Call your function.
xy = sigmaEllipse2D(mu, P);

%Now plot the generated points. You should see an elongated ellipse stretching from the top left corner to the bottom right. 
figure(1);
h1 = plot(xy(1,:), xy(2,:));
%Set the scale of x and y axis to be the same. This should be done if the two variables are in the same domain, e.g. both are measured in meters.
axis equal
hold on
%Also plot a star where the mean is, and make it have the same color as the ellipse.
plot(mu(1), mu(2), '*', 'color', h1.Color);


%% Test approxGaussianTransform

f = @(x)[x(2,:).*cos(x(1,:)); x(2,:).*sin(x(1,:))];
mu_x = [3;-2];
Sigma_x = [2, 0.3; 0.3, 1];
N = 100000;
[mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N);

mu_y
Sigma_y

%% Test jointGaussian

mu_x = 19;
sigma2_x = 5^2;
sigma2_r = 2^2;

[mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r)

%% test posteriorGaussian

mu_x = 3.3575;
sigma2_x = 2.2732;
sigma2_r = 3.4704;
y = 5.4159;
[mu1,sigma1] = posteriorGaussian(mu_x, sigma2_x, y, sigma2_r);
