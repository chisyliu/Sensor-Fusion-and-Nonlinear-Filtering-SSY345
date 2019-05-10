clear all; close all; clc;
cp = fp.getColor(1:10);


%% 1A

sigma_v = 1        *1e-4;
sigma_w = pi/180 ;

name2save = 'S';
savefig = false;

% True track
% Sampling period
T = 0.1;
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(200:400) = -pi/201/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;
% Create true track
for i=2:K+1
    % Simulate
    X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
    % Set turn−rate
    X(5,i) = omega(i);
end


% Prior information
x_0 = [0 0 0 0 0]';
P_0 = diag([10 10 10 5*pi/180 pi/180].^2);
% Sensor positions
sk = [280 -140]';


% measurement noise
R = diag([15 2*pi/180].^2);
% generate measurement sequence
h = @(x) rangeBearingMeasurements(x, sk);
Y = genNonLinearMeasurementSequence(X,h,R);


% Motion model
f = @(x) coordinatedTurnMotion(x,T);
Q = diag([0 0 T*sigma_v^2 0 T*sigma_w^2]);

% [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, 'CKF');
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, Q, h, R, @sigmaPoints, 'CKF');

% calcualte unfiltered position from sensors given angles
Xm = sk + Y(1,:).*[cos(Y(2,:));sin(Y(2,:))];


figure('Color','white','Position',[387  228  990  357]);
subplot(1,2,1)
title('Filter')
plotTurnU( X, xf, Pf, Xm, sk, 'filter', 2)
subplot(1,2,2)
title('Smoother')
plotTurnU( X, xs, Ps, Xm, sk, 'smoother', 4)

figure('Color','white','Position',[428  692  930  207]);
plotTurnUError( T, xf, xs, X )


%% 1B
close all;

Y(1,150) = norm(X(1:2,150)) + 100;


% [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, 'CKF');
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, Q, h, R, @sigmaPoints, 'CKF');

% calcualte unfiltered position from sensors given angles
Xm = sk + Y(1,:).*[cos(Y(2,:));sin(Y(2,:))];

figure('Color','white','Position',[387  228  990  357]);
subplot(1,2,1)
title('Filter')
plotTurnU( X, xf, Pf, Xm, sk, 'filter', 2)
subplot(1,2,2)
title('Smoother')
plotTurnU( X, xs, Ps, Xm, sk, 'smoother', 4)

figure('Color','white','Position',[428  692  930  207]);
plotTurnUError( T, xf, xs, X )


%% 2A


close all;

T = 0.01;

% Motion Parameters
A = 1;
f = @(x) A*x;
Q = 1.5;
% Measurement Parameters
H = 1;
h = @(x) H*x;
R = 2.5;
% Prior
x_0 = 2;
P_0 = 6;

% General parameters
N = 50;

% calculate state and measurement sequences
X = genLinearStateSequence(x_0,P_0,A,Q,N);
Y = genLinearMeasurementSequence(X, H, R);

% filter data
[xf, Pf] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);
% [xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, Q, h, R, @sigmaPoints, 'CKF');
figure
[xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, f, Q, h, R, 1000, true, @plotPostPdf);

plotPostPdf(k, Xk, Wk, xf, Pf, bResample, sigma, ax)

% plot position
figure('Color','white','Position',[199   288  1152   457]);
hold on, grid on;
p3 = plot(0:N, H*X, 'b', 'LineWidth',3, 'DisplayName','true state');
p3.Color = [p3.Color 0.2];
p2 = plot(0:N, H*[x_0 xf], 'b', 'LineWidth',1.5, 'DisplayName','state estimate');
p1 = plot(1:N, Y, '*r', 'DisplayName','measurements');
p4 = plot(0:N, H*[x_0 xf] + 3*sqrt([P_0(1) squeeze(Pf(1,1,:))']), '--b', 'DisplayName','+3-sigma level');
p5 = plot(0:N, H*[x_0 xf] - 3*sqrt([P_0(1) squeeze(Pf(1,1,:))']), '--b', 'DisplayName','-3-sigma level');
xlabel('k - time step');
ylabel('position');
legend([p1 p2 p3 p4 p5],'Location','northwest');
% ylim([-3 11])
% fp.savefig('q2b-pos')



%%


function plotTurnU( X, xf, Pf, Xm, sk, signame, coln)
    cp = fp.getColor(1:10);
    grid on; hold on, axis equal;
    for i=1:5:length(xf)
        ell_xy = sigmaEllipse2D(xf(1:2,i),Pf(1:2,1:2,i),3,50);
        p5 = fill(ell_xy(1,:),ell_xy(2,:), cp(coln,:),'facealpha',.1, 'DisplayName',[signame,' 3-sigma']);   %,'edgecolor','none'
    end

    p1 = plot(X(1,:),X(2,:),   '-', 'Color', cp(1,:), 'LineWidth',3, 'DisplayName','True position sequence');
    p2 = plot(xf(1,:),xf(2,:), '-', 'Color', cp(coln,:), 'LineWidth',3, 'DisplayName',[signame,' position']);
    sc1 = scatter(sk(1), sk(2), 100, 'o', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor', cp(4,:), 'MarkerEdgeColor', cp(4,:),'DisplayName','sensor 1 location');

    axis manual
    p4 = plot(Xm(1,:),Xm(2,:), 'Color', [cp(3,:) 0.3], 'LineWidth',1, 'DisplayName','Measured position');

    xlabel 'pos x', ylabel 'pos y'
    legend([p1 p2 p4 sc1 p5], 'Location','west')
    % if savefig fp.savefig(sprintf('q3_%s',name2save)); end
end

function plotTurnUError( T, xf, xs, X )
    cp = fp.getColor(1:10);
    K = length(X)-1;
    grid on, hold on;
    p1 = plot( (1:K)*T, vecnorm(xf(1:2,:)-X(1:2,2:end), 2, 1), 'Color', cp(2,:) , 'LineWidth',2, 'DisplayName','Filter error');
    p2 = plot( (1:K)*T, vecnorm(xs(1:2,:)-X(1:2,2:end), 2, 1), 'Color', cp(4,:) , 'LineWidth',2, 'DisplayName','Smoother error');
    ylabel('$|p_k - \hat{p_{k|k}}|_2$', 'Interpreter','Latex', 'FontSize',16), xlabel('Time [s]')
    title 'Position error'
    legend([p1 p2])
    % if savefig fp.savefig(sprintf('q3_%s_err',name2save)); end
end

