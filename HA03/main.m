close all; clear all; clc;

cp = fp.getColor(1:10);

%% Approximations of mean and covariance

N=10000;
type = 'EKF';
distribution = 2; % 1 or 2



if distribution == 1
    x_mu = [120 120]';
    x_sigma  = diag([5^2 10^2]);
else
    x_mu = [120 -20]';
    x_sigma  = diag([5^2 10^2]);
end
s1 = [0, 100]';
s2 = [100, 0]';
R = diag([0.1*pi/180 0.1*pi/180]);
h = @(x) dualBearingMeasurement(x,s1,s2);



% A

[y_mu, y_sigma, y, x] = approxGaussianTransform( x_mu, x_sigma, @(x)genNonLinearMeasurementSequence(x,h,R) , N );

if strcmp(type,'UKF') || strcmp(type,'CKF')
    [SP1,W1] = sigmaPoints(x_mu,x_sigma,type)
    hSP1 = h(SP1)
    [ye_mu, ye_sigma]  = estimateSP(h, R, SP1, W1)
else
    [hx, dhx] = h(x_mu);
    % predict using first order Taylor expansion
    ye_mu = hx;
    ye_sigma = dhx * x_sigma * dhx' + R;
end


figure('Color','white','Position',[837  424  603  429]);
grid on; hold on %, axis equal
sc1 = scatter(y(1,:), y(2,:), 20, 'filled', 'MarkerFaceColor', cp(1,:), 'MarkerFaceAlpha',0.1, 'DisplayName','y = h(x) + r');

[ xy ] = sigmaEllipse2D( y_mu, y_sigma, 3, 100 );
p1 = plot(xy(1,:),xy(2,:), 'Color', cp(2,:), 'LineWidth',2, 'DisplayName','Sample 3-sigma ellipse');
sc2 = scatter(y_mu(1), y_mu(2), 100, 'o', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor', cp(2,:), 'MarkerEdgeColor', cp(2,:),'DisplayName','Sample mean');
% sc2 = scatter(y_mu(1), y_mu(2), 100, 'h','filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', cp(2,:), 'MarkerEdgeColor', cp(2,:),'DisplayName','Sample mean');

if strcmp(type,'UKF') || strcmp(type,'CKF')
    sc3 = scatter(hSP1(1,:), hSP1(2,:), 100, 'h','filled', 'MarkerFaceAlpha',1, 'MarkerFaceColor', cp(4,:), 'MarkerEdgeColor', cp(4,:),'DisplayName','Sigma points');
end

sc4 = scatter(ye_mu(1,:), ye_mu(2,:), 100, 'o','filled', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor', cp(4,:), 'MarkerEdgeColor', cp(4,:),'DisplayName','Approximated mean');
% sc4 = scatter(ye_mu(1,:), ye_mu(2,:), 100, 'h','filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', cp(4,:), 'MarkerEdgeColor', cp(4,:),'DisplayName','Approximated mean');
[ xy ] = sigmaEllipse2D( ye_mu, ye_sigma, 3, 100 );
p2 = plot(xy(1,:),xy(2,:), '--', 'Color', cp(4,:), 'LineWidth',3,'DisplayName','Approximated 3-sigma ellipse');

title(sprintf('Type: %s, x~p%d(x)',type,distribution))

legend('Location','southeast')

fp.savefig(sprintf('q1_t_%s_d_%d',type,distribution))



%% Non-linear Kalman filtering





%% Help functions

function [x, P] = estimateSP(f, R, SP, W)
    n = size(SP,1);
    x = zeros(n,1);
    for i=1:numel(W)
        x = x + f(SP(:,i)) * W(i);
    end   
    P = R; %zeros(n,n);
    for i=1:numel(W)
        P = P + (f(SP(:,i))-x)*(f(SP(:,i))-x).' * W(i);
    end
end


