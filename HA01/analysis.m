clearvars; close all; clc;

%% Transformation of Gaussian random variables
% (a)

close all;

mu_x = [10;0];
Sigma_x = diag([0.2 8]);

A = [1 1; 1 -1];
b = [0;0];
[mu_y, Sigma_y] = affineGaussianTransform(mu_x,Sigma_x, A, b);
xy = sigmaEllipse2D(mu_y,Sigma_y,3,100);


f = @(x) A*x+b;
N = 1000;
[mu_ya, Sigma_ya, y_s] = approxGaussianTransform(mu_x,Sigma_x,f,N);
xya = sigmaEllipse2D(mu_ya,Sigma_ya,3,100);



figure('Color','white','Position',[670  384  772  580]); hold on, grid on;
cp = fp.getColor(1:2);

scp = scatter(y_s(1,:), y_s(2,:), 10,'filled', 'MarkerFaceColor',cp(1,:));
scp.MarkerFaceAlpha = 0.5;

h1a = plot(xya(1,:), xya(2,:), 'LineWidth',3, 'Color', cp(1,:));
s1a = scatter(mu_ya(1), mu_ya(2), 100, 'h','filled', 'MarkerFaceColor', cp(1,:), 'MarkerEdgeColor', cp(1,:));

h1 = plot(xy(1,:), xy(2,:),'--', 'LineWidth',3, 'Color', cp(2,:));
sc1= scatter(mu_y(1), mu_y(2), 100, 'h','filled', 'MarkerFaceColor', cp(2,:), 'MarkerEdgeColor', cp(2,:));
sc1.MarkerFaceAlpha = 1;

xlim([0,20])
ylim([0,20])

xlabel 'y(1)', ylabel 'y(2)'
title({'Approximated and analytical computation of transformed density p(y)=h(x)',['N=',num2str(N)]})
legend({'Approx transformed samples', 'Approx. 3\sigma ellipse', 'Approx. mean', ...
        'Analytic  3\sigma ellipse','Analytic mean'})
% fp.savefig(sprintf('q1-N-%s',num2str(N)))



%% Transformation of Gaussian random variables 
%  (b)

close all;

mu_x = [10;0];
Sigma_x = diag([0.2 8]);

f = @(x) [ (x(1,:).^2+x(2,:).^2).^0.5 ; atan2(x(2,:),x(1,:))];
N = 10;
[mu_ya, Sigma_ya, y_s] = approxGaussianTransform(mu_x,Sigma_x,f,N);
xya = sigmaEllipse2D(mu_ya,Sigma_ya,3,100);


figure('Color','white','Position',[665  462  670  497]);
hold on, grid on;
cp = fp.getColor(1:2);

scp = scatter(y_s(1,:), y_s(2,:), 10,'filled', 'MarkerFaceColor',cp(1,:));
scp.MarkerFaceAlpha = 0.3;

h1a = plot(xya(1,:), xya(2,:), 'LineWidth',3, 'Color', cp(1,:));
s1a = scatter(mu_ya(1), mu_ya(2), 100, 'h','filled', 'MarkerFaceColor', cp(1,:), 'MarkerEdgeColor', cp(1,:));

xlim([8,13])
ylim([-1,1])

xlabel 'y(1)', ylabel 'y(2)'
title({'Approximated and analytical computation of transformed density p(y)=h(x)',['N=',num2str(N)]})
legend({'Approx transformed samples', 'Approx. 3\sigma ellipse', 'Approx. mean'})
% fp.savefig(sprintf('q1b-N-%s',num2str(N)))



%% Snow depth in Norway
% (a)

close all;

mu_x = 1.1;
sigma2_x = 0.5^2;
sigma2_r = 0.2^2;
[mu_xy, sigma_xy] = jointGaussian(mu_x, sigma2_x, sigma2_r);
xy = sigmaEllipse2D(mu_xy,sigma_xy,3,100);


figure('Color','white','Position',[665  462  670  497]);
hold on, grid on, axis equal;
cp = fp.getColor(1:2);

% scp = scatter(y_s(1,:), y_s(2,:), 10,'filled', 'MarkerFaceColor',cp(1,:));
% scp.MarkerFaceAlpha = 0.3;

h1a = plot(xy(1,:), xy(2,:), 'LineWidth',3, 'Color', cp(1,:));
s1a = scatter(mu_xy(1), mu_xy(2), 100, 'h','filled', 'MarkerFaceColor', cp(1,:), 'MarkerEdgeColor', cp(1,:));

xlim([-1,3])
ylim([-1,3])

xlabel 'x', ylabel 'y'% scp = scatter(y_s(1,:), y_s(2,:), 10,'filled', 'MarkerFaceColor',cp(1,:));
% scp.MarkerFaceAlpha = 0.3;

title({'Joint distribution p(x,y) for Anna'})
legend({'3\sigma ellipse', 'Mean'})
fp.savefig(sprintf('q2a-Anna'))



mu_x = 1;
sigma2_x = 0.5^2;
sigma2_r = 1^2;
[mu_xy, sigma_xy] = jointGaussian(mu_x, sigma2_x, sigma2_r);
xy = sigmaEllipse2D(mu_xy,sigma_xy,3,100);


figure('Color','white','Position',[359  250  384  497]);
hold on, grid on, axis equal;
cp = fp.getColor(1:2);

% scp = scatter(y_s(1,:), y_s(2,:), 10,'filled', 'MarkerFaceColor',cp(1,:));
% scp.MarkerFaceAlpha = 0.3;

h1a = plot(xy(1,:), xy(2,:), 'LineWidth',3, 'Color', cp(1,:));
s1a = scatter(mu_xy(1), mu_xy(2), 100, 'h','filled', 'MarkerFaceColor', cp(1,:), 'MarkerEdgeColor', cp(1,:));

xlim([-1 3])
ylim([-3,5])

xlabel 'x', ylabel 'y'% scp = scatter(y_s(1,:), y_s(2,:), 10,'filled', 'MarkerFaceColor',cp(1,:));
% scp.MarkerFaceAlpha = 0.3;

title({'Joint distribution p(x,y) for Else'})
legend({'3\sigma ellipse', 'Mean'})
fp.savefig(sprintf('q2a-Else'))



%% Snow depth in Norway
% (b)


%% Snow depth in Norway
% (c)







%% MMSE and MAP estimates for Gaussian mixture posteriors


